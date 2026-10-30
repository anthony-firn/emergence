"""Browser QA for the cross-domain workspace: command line, tour, links,
sensitivities, every domain studio, 3D render and verified exports."""
import argparse
import json
import zipfile
from pathlib import Path
from playwright.sync_api import sync_playwright
from test_cad_browser import POLICY

DOCS = "() => JSON.stringify(window.emergence.documents())"


def check(url, screenshots=None):
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch(args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
        context = browser.new_context(viewport={"width": 1440, "height": 1000}, color_scheme="dark", accept_downloads=True)

        def policy(route):
            response = route.fetch()
            route.fulfill(response=response, headers={**response.headers, "content-security-policy": POLICY})

        context.route("**/play/emergence/cad/", policy)
        page = context.new_page()
        page.on("pageerror", lambda error: errors.append(str(error)))
        solved = "() => window.emergence.solution && window.emergence.solution.outputs['flight.total_mass']"
        page.goto(url, wait_until="networkidle")
        page.wait_for_function(solved, timeout=60000)
        assert page.locator("#system-studio").is_visible()
        value = lambda port: page.evaluate(f"() => window.emergence.solution.outputs['{port}'].value")

        def settle():
            page.wait_for_function("() => !document.body.classList.contains('solving')", timeout=60000)

        mass, current = value("flight.total_mass"), value("flight.esc_current")
        assert page.evaluate("() => window.emergence.solution.converged")
        # The command line edits a domain and the linked solve carries it across.
        page.locator("#command-input").fill("flight.mass gimbal 1.4")
        page.keyboard.press("Enter")
        page.wait_for_function("() => document.querySelector('#command-output').textContent.startsWith('✓')")
        page.wait_for_function(f"() => window.emergence.solution.outputs['flight.total_mass'].value > {mass + .8}", timeout=60000)
        assert value("flight.esc_current") > current
        page.locator("#command-input").fill("undo")
        page.keyboard.press("Enter")
        page.wait_for_function(f"() => Math.abs(window.emergence.solution.outputs['flight.total_mass'].value - {mass}) < 1e-9", timeout=60000)
        page.locator("#command-input").fill("pcb.plac")
        assert page.locator("#command-suggestions li").first.inner_text().startswith("pcb.place")
        page.locator("#command-input").fill("")
        print("Command line, linked propagation and undo passed.")

        before = page.evaluate(DOCS)
        page.locator("#start-tour").click()
        steps = 0
        while True:
            label = page.locator("#tour-next").inner_text()
            page.locator("#tour-next").click()
            if label == "Run step":
                page.wait_for_function("() => !document.querySelector('#tour-next').disabled && !document.querySelector('#tour-result').textContent.startsWith('Running')", timeout=90000)
                result = page.locator("#tour-result").inner_text()
                assert "→" in result, result
                steps += 1
                continue
            if label == "Finish":
                break
        assert steps >= 6 and page.locator("#tour").is_hidden()
        settle()
        assert page.evaluate(DOCS) == before, "The tour restores the design exactly"
        print(f"Guided tour ran {steps} live cross-domain commands and restored the design.")

        # Links: disabling one removes the drive; the System graph shows every domain.
        assert page.locator(".link-graph .domain-node").count() >= 8
        page.locator(".link-table tr", has_text="Design-point ESC current").locator("input[type=checkbox]").uncheck()
        page.wait_for_function("() => !window.emergence.solution.applied['mech.current']", timeout=60000)
        page.locator("#command-undo").click()
        page.wait_for_function("() => window.emergence.solution.applied['mech.current']", timeout=60000)
        page.locator("#study-output").select_option("mech.q1_tj")
        page.locator("#study-run").click()
        page.wait_for_function("() => document.querySelectorAll('.tornado > div').length === 8", timeout=240000)
        print("Link toggling, undo and the coupled sensitivity study passed.")

        page.locator("#studio-tab-optics").click()
        assert page.locator(".spot").count() == 3
        z = page.evaluate("() => window.emergence.documents().optics.target.z")
        page.locator("#optics-focus").click()
        page.wait_for_function(f"() => window.emergence.documents().optics.target.z !== {z}")
        page.locator("#command-undo").click()
        page.locator("#studio-tab-parts").click()
        page.wait_for_function("() => document.querySelector('#parts-mesh-info').dataset.mesh === 'ready'", timeout=90000)
        assert "watertight" in page.locator("#parts-mesh-info").inner_text()
        page.locator("#studio-tab-render").click()
        page.wait_for_function("() => / objects · .* triangles · from /.test(document.querySelector('#render-status').textContent)", timeout=120000)
        with page.expect_download() as png:
            page.locator("#render-png").click()
        assert Path(png.value.path()).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        if screenshots:
            screenshots.mkdir(parents=True, exist_ok=True)
            page.locator("#render-stage, .render-stage").first.screenshot(path=str(screenshots / "emergence-cad-render.png"))
        print("Optics best focus, kernel part meshes and the whole-system render passed.")

        page.locator("#studio-tab-pcb").click()
        with page.expect_download(timeout=120000) as fab:
            page.locator("#pcb-fab").click()
        archive = zipfile.ZipFile(fab.value.path())
        assert {"emergence-F_Cu.gbr", "emergence-B_Cu.gbr", "emergence.drl", "emergence.kicad_pcb"} <= set(archive.namelist())
        assert archive.read("emergence.drl").startswith(b"M48")
        page.locator("#studio-tab-silicon").click()
        with page.expect_download(timeout=60000) as gds:
            page.locator("#silicon-gds").click()
        assert Path(gds.value.path()).read_bytes()[:4] == b"\x00\x06\x00\x02"
        assert "read back exactly" in page.locator("#silicon-status").inner_text()
        page.locator("#litho-ilt").click()
        page.wait_for_function("() => document.querySelector('#litho-status').textContent.startsWith('ILT')", timeout=60000)
        page.wait_for_selector("#silicon-results .critical-path")
        print("Fabrication package, verified GDSII, inverse lithography and the chip flow passed.")

        # One save carries every domain's document into the revision graph.
        page.locator("#command-input").fill("architecture.wall charge west exterior=true window=4")
        page.keyboard.press("Enter")
        page.wait_for_function("() => !document.querySelector('#global-save').disabled")
        page.locator("#global-save").click()
        page.locator("#save-message").fill("Hangar: west glazing")
        page.locator('#save-form button[type="submit"]').click()
        page.wait_for_function("() => document.querySelector('#revision-count').textContent === '3'")
        with page.expect_download() as exported:
            page.locator("#export").click()
        bundle = json.loads(Path(exported.value.path()).read_text())
        assert bundle["integrity"]["algorithm"] == "SHA-256" and set(bundle["commits"][-1]["studios"]) >= {"flight", "pcb", "silicon", "optics", "fluids", "architecture", "parts", "system"}
        if screenshots:
            for studio in ("system", "flight", "parts", "optics", "fluids", "architecture"):
                page.locator(f"#studio-tab-{studio}").click()
                page.evaluate("scrollTo({top:0,behavior:'instant'}); document.querySelector('#toast').hidden=true")
                page.screenshot(path=str(screenshots / f"emergence-cad-{studio}.png"))
        browser.close()
    assert not errors, errors
    print("Global revision of every domain with a verified, sealed export passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url", nargs="?", default="http://127.0.0.1:18762/play/emergence/cad/")
    parser.add_argument("--screenshots", type=Path)
    args = parser.parse_args()
    check(args.url, args.screenshots)
