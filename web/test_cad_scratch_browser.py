"""Browser QA: a person starts from nothing. Every edit below goes through a
visible control, never the command line or the page's API, which are only read
to check the result. A new design is built across every domain from a blank
project, then the drone's build script is replayed into another blank project
and must give back the shipped design."""
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

        def fill(selector, value):
            page.locator(selector).fill(str(value))
            page.locator(selector).press("Tab")

        def tab(name):
            page.locator(f"#studio-tab-{name}").click()

        def docs():
            return page.evaluate("() => window.emergence.documents()")

        def solved(expression, timeout=60000):
            page.wait_for_function(expression, timeout=timeout)

        def add_domain(name):
            page.locator("#studio-add").click()
            page.locator(f'#domain-choices [data-domain="{name}"]').get_by_role("button", name="Add empty").click()
            page.wait_for_function(f"() => !document.querySelector('#{name}-studio').hidden")

        def form(prefix, values, submit=True):
            for key, value in values.items():
                field = page.locator(f"#{prefix}-{key}")
                if field.evaluate("e => e.tagName") == "SELECT":
                    field.select_option(str(value))
                else:
                    field.fill(str(value))
            if submit:
                page.locator(f"#{prefix}-submit").click()

        # Nothing: no domains, no links, no documents to speak of.
        page.goto(url + "?new=blank#system", wait_until="networkidle")
        page.wait_for_selector(".empty-project", timeout=60000)
        page.locator("#project-edit").click()
        page.locator("#project-title").fill("Handheld inspection lamp")
        page.locator("#project-description").fill("A battery lamp with a focusing lens, its driver board and a cooled housing.")
        page.locator('#project-form button[type="submit"]').click()
        page.wait_for_function("() => document.querySelector('#project-name').textContent === 'Handheld inspection lamp'")

        # A housing, built feature by feature from a disc.
        page.locator("#system-add-parts").click()
        tab("parts")
        assert "No parts" in page.locator("#parts-studio").inner_text() or page.locator("#parts-studio .studio-tree-item").count() == 0
        page.locator("#parts-template").select_option("disc")
        page.locator("#parts-add").click()
        page.wait_for_selector("#parts-define-name")
        form("parts-define", {"name": "bore", "value": "diameter/3"})
        solved("() => window.emergence.documents().parts.parts[0].model.params.bore === 'diameter/3'")
        form("parts-solid", {"how": "cut", "shape": "cylinder", "radius": "bore/2", "height": 40, "z": 6})
        solved("() => window.emergence.documents().parts.parts[0].model.shape.op === 'difference'")
        form("parts-solid", {"how": "add", "shape": "extrude", "points": "[[28,-6],[44,-6],[44,6],[28,6]]", "height": 8, "z": 4})
        solved("() => JSON.stringify(window.emergence.documents().parts.parts[0].model.shape).includes('extrude')")
        form("parts-modify", {"op": "polar", "count": 3})
        solved("() => window.emergence.documents().parts.parts[0].model.shape.op === 'polar_array'")
        page.wait_for_function("() => document.querySelector('#parts-mesh-info').dataset.mesh === 'ready'", timeout=90000)
        assert "NOT watertight" not in page.locator("#parts-mesh-info").inner_text()
        assert page.locator(".feature-outline li").count() >= 4
        fill("#parts-name", "Housing")
        page.locator("#parts-pattern").select_option("stack")
        page.wait_for_selector("#parts-count")
        fill("#parts-count", 2)
        fill("#parts-stack", 14)
        solved("() => { const p = window.emergence.documents().parts.parts[0]; return p.pattern === 'stack' && p.count === 2 && p.stack === 14 && p.name === 'Housing'; }")
        print("A housing was built from a disc with a parametric bore, a sketch, a polar pattern and a stack.")

        # A driver board on a footprint that is not built in.
        add_domain("pcb")
        page.locator(".footprint-tools summary").click()
        form("pcb-footprint", {"name": "soic14", "pins": 14, "pitch": 1.27, "span": 5.4, "pad": .8, "drill": .4})
        page.wait_for_function("() => [...document.querySelectorAll('#pcb-type option')].some(o => o.value === 'soic14')")
        page.locator("#pcb-type").select_option("soic14")
        page.locator("#pcb-add").click()
        page.wait_for_selector("#pcb-net-14")
        fill("#pcb-name", "LED driver")
        fill("#pcb-net-1", "VIN")
        for pad in range(2, 15):
            fill(f"#pcb-net-{pad}", "NC")
        page.locator("#pcb-type").select_option("connector")
        page.locator("#pcb-add").click()
        solved("() => window.emergence.documents().pcb.components.length === 2")
        page.wait_for_selector("#pcb-net-2")
        fill("#pcb-net-1", "VIN")
        fill("#pcb-net-2", "NC")
        fill("#pcb-x", 12)
        solved("() => window.emergence.documents().pcb.components[1].nets[0] === 'VIN'")
        page.locator("#pcb-route-all").click()
        solved("() => window.emergence.documents().pcb.routes.length === 1")
        page.locator("#pcb-connectivity").click()
        fill("#pcb-current-VIN", .5)
        solved("() => window.emergence.documents().pcb.currents.VIN === 0.5")
        solved("() => { const s = window.emergence.solution; return s?.results.pcb?.details.board.routes.length === 1 && s.results.pcb.details.board.currents.VIN === 0.5 && s.checks.some(c => c.domain === 'pcb' && c.id === 'geometry' && c.status === 'met'); }")
        print("A custom footprint was made, placed, netted, routed and given a current; its geometry checks pass.")

        # A lens, chosen by hand.
        add_domain("optics")
        page.get_by_role("button", name="＋ Lens").click()
        page.wait_for_selector("#optics-L1-r1")
        fill("#optics-L1-r1", 40)
        fill("#optics-L1-r2", -40)
        fill("#optics-L1-thickness", 6)
        fill("#optics-L1-name", "Front lens")
        fill("#optics-fields", "0,3")
        fill("#optics-power", 2)
        page.locator("#optics-focus").click()
        page.wait_for_function("() => { const d = window.emergence.documents().optics; return d.elements[0].name === 'Front lens' && d.fields.length === 2 && Math.abs(d.target.z - 100) > 1; }")

        # Cooling air through fins, in a room, with one formula of our own.
        add_domain("fluids")
        page.get_by_role("button", name="＋ Duct").click()
        page.locator("#fluids-add-channels").click()
        solved("() => window.emergence.documents().fluids.elements.length === 2")
        add_domain("architecture")
        page.wait_for_function("() => document.querySelector('#architecture-studio').innerText.includes('No rooms yet')")
        page.get_by_role("button", name="＋ Add room").click()
        page.wait_for_selector("#arch-name")
        fill("#arch-name", "Test bench")
        fill("#arch-site-equipment", 40)
        solved("() => window.emergence.documents().architecture.rooms[0].name === 'Test bench' && window.emergence.documents().architecture.equipmentPower === 40")
        add_domain("equations")
        form("equations-input", {"id": "board_power", "value": 1, "unit": "W", "label": "Board power"})
        page.wait_for_selector("#equations-value-board_power")
        form("equations-input", {"id": "battery_wh", "value": 7.4, "label": "Battery energy"})
        page.wait_for_selector("#equations-value-battery_wh")
        form("equations-output", {"id": "hours", "expr": "battery_wh / (board_power + 0.4)", "unit": "h", "label": "Runtime"})
        solved("() => window.emergence.solution?.outputs['equations.hours']?.value > 5")
        print("Optics, fluids, a room and a formula were all started empty.")

        # Tie it together: a link, a requirement, a headline, a study, an optimisation, a tour, a credit.
        tab("system")
        page.wait_for_selector(".link-form")
        page.locator(".link-form select").first.select_option("pcb.dissipation")
        page.locator(".link-form select").nth(1).select_option("equations.board_power")
        page.locator("#system-add-link").click()
        solved("() => window.emergence.solution?.applied['equations.board_power']")
        page.wait_for_selector("#requirement-port")
        form("requirement", {"port": "equations.hours", "min": 2, "label": "Runs at least two hours"})
        solved("() => window.emergence.solution?.checks.some(c => c.domain === 'project' && c.status === 'met')")
        page.select_option("#system-pin", "equations.hours")
        page.locator("#system-pin-add").click()
        page.wait_for_function("() => window.emergence.documents().project.headline.length === 1")
        page.wait_for_selector("#study-domain")
        form("study", {"domain": "equations", "path": "inputs.battery_wh.value", "label": "Battery energy"})
        page.wait_for_function("() => window.emergence.documents().project.studies.length === 1")
        page.locator("#study-run").click()
        page.wait_for_selector(".tornado", timeout=60000)
        form("optimize", {"title": "Smallest battery", "objective": "equations.hours", "sense": "min"})
        page.wait_for_selector("#optimize-variable-domain")
        page.locator("#optimize-variable-domain").select_option("equations")
        form("optimize-variable", {"path": "inputs.battery_wh.value", "min": 1, "max": 20, "step": .5, "label": "Battery energy"})
        page.wait_for_function("() => window.emergence.documents().project.optimize.variables.length === 1")
        page.wait_for_selector("#optimize-constraint-port")
        form("optimize-constraint", {"port": "equations.hours", "min": 2, "label": "At least two hours"})
        page.wait_for_function("() => window.emergence.documents().project.optimize.constraints.length === 1")
        page.locator("#optimize-run").click()
        page.wait_for_selector(".optimum", timeout=120000)
        assert page.locator(".optimum").inner_text().startswith("Best after"), page.locator(".optimum").inner_text()
        page.locator("#optimize-apply").click()
        page.wait_for_function("() => window.emergence.documents().equations.inputs.find(x => x.id === 'battery_wh').value < 7.4")
        page.wait_for_selector("#tour-step-form-title")
        form("tour-step-form", {"title": "Meet the lamp", "text": "One file, seven domains."})
        page.wait_for_function("() => window.emergence.documents().project.tour.length === 1")
        page.wait_for_selector("#start-tour:not([hidden])")
        page.locator("details.authoring summary", has_text="Credit a design").click()
        form("credit", {"title": "Plano-convex collimator note", "url": "https://example.org/lens", "licence": "CC0"})
        page.wait_for_function("() => window.emergence.documents().project.credits.length === 1")
        print("A link, requirement, headline, study, optimisation, tour step and credit were authored.")

        # The 3D scene shows every domain, placed by hand.
        tab("render")
        page.wait_for_function("() => / triangles · from /.test(document.querySelector('#render-status')?.textContent || '')", timeout=120000)
        fill("#render-pcb-z", 10)
        page.wait_for_function("() => window.emergence.documents().project.scene.placements.pcb?.at[2] === 10")
        status = page.locator("#render-status")
        page.wait_for_function("() => ['Parts', 'PCB', 'Optics'].every(d => document.querySelector('#render-status').textContent.includes(d))", timeout=120000)
        print("Render:", status.inner_text()[:90])

        # Every authoring panel fits a phone.
        for width in (320, 390, 768):
            page.set_viewport_size({"width": width, "height": 844})
            for studio in ("system", "parts", "pcb", "optics", "fluids", "architecture", "equations", "render"):
                tab(studio)
                page.wait_for_timeout(120)
                overflow = page.evaluate("document.documentElement.scrollWidth - innerWidth")
                assert overflow <= 1, (width, studio, overflow)
        page.set_viewport_size({"width": 1440, "height": 1000})
        print("Every authoring panel fits 320, 390 and 768 px.")

        # It survives a reload as a saved revision.
        page.wait_for_function("() => !document.querySelector('#global-save').disabled")
        page.locator("#global-save").click()
        page.locator("#save-message").fill("First complete lamp")
        page.locator('#save-form button[type="submit"]').click()
        page.wait_for_function("() => document.querySelector('#global-save').disabled")
        before = docs()
        page.reload(wait_until="networkidle")
        page.wait_for_function("() => window.emergence?.solution", timeout=120000)
        after = docs()
        for key in ("project", "system", "parts", "pcb", "optics", "fluids", "architecture", "equations"):
            assert after[key] == before[key], key
        assert page.locator("#project-name").inner_text() == "Handheld inspection lamp"
        if screenshots:
            screenshots.mkdir(parents=True, exist_ok=True)
            tab("system")
            page.wait_for_selector(".system-hero")
            page.evaluate("document.querySelector('#command-dock').style.display='none'; scrollTo({top:0,behavior:'instant'})")
            page.screenshot(path=str(screenshots / "emergence-cad-from-scratch.png"))
            page.evaluate("document.querySelector('#command-dock').style.display=''")
        print("The lamp was saved as a revision and reloaded identically.")

        # The drone is reachable the same way: its own build script, replayed into a blank project.
        page.goto(url + "?example=drone#system", wait_until="networkidle")
        page.wait_for_function("() => window.emergence?.solution && document.querySelector('#system-build-script')", timeout=120000)
        shipped = docs()
        with page.expect_download() as dl:
            page.locator("#system-build-script").click()
        script = Path(dl.value.path()).read_bytes()
        assert len(json.loads(script)) > 100
        page.goto(url + "?new=blank#system", wait_until="networkidle")
        page.wait_for_selector("#system-replay", state="attached", timeout=60000)
        page.locator("#system-replay").set_input_files({"name": "drone-build.json", "mimeType": "application/json", "buffer": script})
        page.wait_for_function("() => document.querySelector('#project-name').textContent === 'Inspection quadcopter'", timeout=60000)
        page.wait_for_function("() => window.emergence.documents().project.tour.length === 13")
        rebuilt = docs()
        for key in ("project", "system", "flight", "mech", "parts", "structures", "power", "pcb", "silicon", "litho", "optics", "fluids", "architecture", "controls", "testing"):
            assert rebuilt[key] == shipped[key], key
        print("The drone's build script replayed into a blank project and gave back the shipped design.")
        browser.close()
    assert not errors, errors
    print("From nothing, through visible controls alone, passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url", nargs="?", default="http://127.0.0.1:18762/play/emergence/cad/")
    parser.add_argument("--screenshots", type=Path)
    args = parser.parse_args()
    check(args.url, args.screenshots)
