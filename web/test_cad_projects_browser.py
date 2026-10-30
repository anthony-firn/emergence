"""Browser QA for projects: examples as project files, the shelter's own
optimisation, project files round-tripping into new projects, a blank project
grown by the UI, dormant links, and histories that never cross projects."""
import argparse
import json
from pathlib import Path
from playwright.sync_api import sync_playwright
from test_cad_browser import POLICY


def check(url, screenshots=None):
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        ctx = browser.new_context(viewport={"width": 1440, "height": 1000}, color_scheme="dark", accept_downloads=True)

        def policy(route):
            response = route.fetch()
            route.fulfill(response=response, headers={**response.headers, "content-security-policy": POLICY})
        ctx.route("**/play/emergence/cad/**", lambda route: policy(route) if route.request.resource_type == "document" else route.continue_())
        ctx.on("page", lambda page: page.on("pageerror", lambda error: errors.append(str(error))))
        page = ctx.new_page()
        tabs = lambda: page.evaluate("() => [...document.querySelectorAll('[data-studio]')].filter(b => !b.hidden).map(b => b.dataset.studio)")
        solved = lambda: page.wait_for_function("() => window.emergence?.solution", timeout=120000)

        # A first visit opens the drone example, one project file among others.
        page.goto(url + "#system", wait_until="networkidle"); solved()
        assert page.locator("#project-name").inner_text() == "Inspection quadcopter" and "example-drone" in page.url
        assert set(tabs()) >= {"flight", "mechanical", "parts", "pcb", "silicon", "optics", "fluids", "architecture"}
        drone = page.evaluate("() => window.emergence.file()")
        assert drone["format"] == "emergence-project/1" and drone["documents"]["project"]["domains"][0] == "flight"

        # The shelter is a different file: no flight, no optics, its own headline and optimisation.
        page.locator("#projects-open").click()
        page.locator("#project-examples .project-card", has_text="Edge compute shelter").locator("button").click()
        page.wait_for_function("() => document.querySelector('#project-name')?.textContent === 'Edge compute shelter'", timeout=60000); solved()
        assert "flight" not in tabs() and "optics" not in tabs() and "architecture" in tabs()
        page.wait_for_selector("#optimize-run")
        assert page.locator(".hero-stat span").first.inner_text() == "Shelter air"
        page.locator("#optimize-run").click()
        page.wait_for_selector(".optimum", timeout=600000)
        result = page.locator(".optimum").inner_text()
        nodes = int(result.split("Nodes ")[1].split(" ")[0].rstrip(","))
        assert result.startswith("Best after") and 100 < nodes < 160, result
        page.locator("#optimize-apply").click()
        page.wait_for_function("n => window.emergence.documents().architecture.units === n", arg=nodes)
        page.locator("#command-undo").click()
        page.wait_for_function("() => window.emergence.documents().architecture.units === 24")
        print(f"The shelter's own optimisation fills it with {nodes} nodes; apply and undo passed.")
        if screenshots:
            screenshots.mkdir(parents=True, exist_ok=True)
            page.evaluate("scrollTo({top:0,behavior:'instant'})"); page.screenshot(path=str(screenshots / "emergence-cad-shelter.png"))

        # A downloaded project file opens as a new project with its own history.
        page.locator("#projects-open").click()
        with page.expect_download() as dl:
            page.locator("#project-download").click()
        file = json.loads(Path(dl.value.path()).read_text())
        assert dl.value.suggested_filename == "edge-compute-shelter.emergence.json" and file["documents"]["project"]["name"] == "Edge compute shelter"
        shelter_history = page.evaluate("() => localStorage.getItem('emergence-cad-project-v1:example-edge-shelter')")
        file["documents"]["project"]["name"] = "Shelter, second site"
        file["documents"]["architecture"]["site"]["outdoor"] = 38
        page.locator("#project-file").set_input_files({"name": "second.emergence.json", "mimeType": "application/json", "buffer": json.dumps(file).encode()})
        page.wait_for_function("() => document.querySelector('#project-name')?.textContent === 'Shelter, second site'", timeout=60000); solved()
        assert page.locator("#revision-count").inner_text() == "2" and page.evaluate("() => window.emergence.documents().architecture.site.outdoor") == 38
        # Another project's history is refused rather than merged.
        page.locator("#import").set_input_files({"name": "other.json", "mimeType": "application/json", "buffer": shelter_history.encode()})
        page.wait_for_function("() => document.querySelector('#toast').textContent.includes('belongs to the project')")
        assert page.locator("#revision-count").inner_text() == "2"
        print("Project files round-trip into new projects; histories never cross projects.")

        # A blank project grows by the UI; removing a domain leaves its links dormant.
        page.locator("#projects-open").click(); page.locator("#project-new").click()
        page.wait_for_selector(".empty-project", timeout=60000)
        assert tabs() == ["system", "render", "history"]
        page.locator("#system-add-optics").click(); page.wait_for_function("() => !document.querySelector('#studio-tab-optics').hidden")
        page.locator("#studio-add").click(); page.locator("#domain-choices .project-card", has_text="Fluids").locator("button").click()
        page.wait_for_function("() => !document.querySelector('#studio-tab-fluids').hidden"); page.locator('[data-close="domain-dialog"]').click()
        page.locator("#command-input").fill("system.link optics.power_on_target fluids.heat_load"); page.keyboard.press("Enter")
        page.locator("#command-input").fill("optics.set power=3"); page.keyboard.press("Enter")
        page.wait_for_function("() => window.emergence.solution?.applied['fluids.heat_load']?.value > 0", timeout=60000)
        page.locator("#command-input").fill("project.remove_domain optics"); page.keyboard.press("Enter")
        page.wait_for_function("() => document.querySelector('#studio-tab-optics').hidden && window.emergence.solution?.dormant?.length === 1", timeout=60000)
        page.locator("#command-input").fill("optics.focus"); page.keyboard.press("Enter")
        page.wait_for_function("() => document.querySelector('#command-output').textContent.includes('not in this project')")
        page.locator("#command-undo").click()
        page.wait_for_function("() => !document.querySelector('#studio-tab-optics').hidden")
        page.locator("#studio-tab-render").click()
        page.wait_for_function("() => / triangles · from .*Optics/.test(document.querySelector('#render-status').textContent)", timeout=120000)
        print("A blank project grew by the UI; dormant links, refused domains and undo passed.")

        # Every project is listed; the drone's own history was never touched.
        page.locator("#projects-open").click()
        names = page.locator("#project-list strong").all_inner_texts()
        assert {"Untitled project", "Shelter, second site", "Edge compute shelter", "Inspection quadcopter"} <= set(names), names
        assert page.evaluate("() => JSON.parse(localStorage.getItem('emergence-cad-project-v1:example-drone')).name") == "Inspection quadcopter"
        for width in (320, 390):
            page.set_viewport_size({"width": width, "height": 844})
            assert page.evaluate("document.querySelector('#projects-dialog').scrollWidth <= innerWidth"), width
        browser.close()
    assert not errors, errors
    print("Projects on this device passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url", nargs="?", default="http://127.0.0.1:18762/play/emergence/cad/")
    parser.add_argument("--screenshots", type=Path)
    args = parser.parse_args()
    check(args.url, args.screenshots)
