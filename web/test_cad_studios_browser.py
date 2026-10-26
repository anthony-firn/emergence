"""Browser QA for mechanical profiles, PCB layout and silicon floorplanning."""
import argparse
import json
from pathlib import Path
from playwright.sync_api import sync_playwright
from test_cad_browser import POLICY


def check(url, screenshots=None):
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        context = browser.new_context(viewport={"width": 1360, "height": 1000},
                                      color_scheme="dark", accept_downloads=True)

        def policy(route):
            response = route.fetch()
            route.fulfill(response=response, headers={**response.headers, "content-security-policy": POLICY})

        context.route("**/play/emergence/cad/", policy)
        page = context.new_page()
        page.on("pageerror", lambda error: errors.append(str(error)))
        def capture(kind):
            page.evaluate("scrollTo(0, 0); document.querySelector('#toast').hidden=true")
            page.screenshot(path=str(screenshots / ("emergence-cad-" + kind + ".png")))

        page.goto(url, wait_until="networkidle")
        page.wait_for_selector(".requirement[data-status=undecided]")
        page.locator("#sketch-part").click()
        page.locator("#part-name").fill("Sketch plate")
        page.locator("#sketch-width").fill("20")
        page.locator("#sketch-depth").fill("10")
        page.locator("#sketch-height").fill("5")
        assert "1000.00 mm³" in page.locator("#sketch-preview").inner_text()
        page.locator('#part-form button[type="submit"]').click()
        assert page.locator("#part-width").input_value() == "20"
        assert page.locator("#part-depth").input_value() == "10"
        assert page.locator("#part-height").input_value() == "5"
        assert "rectangle on XY" in page.locator(".feature-detail").inner_text()
        page.get_by_role("button", name="Duplicate extrusion").click()
        assert page.locator("#part-properties > input").input_value() == "Sketch plate copy"
        page.locator("#undo").click()
        assert page.locator(".feature-detail").count() == 1
        for camera in ("top", "front", "right", "iso"):
            page.locator(f'[data-camera="{camera}"]').click()
        page.locator("#discard").click()
        if screenshots:
            screenshots.mkdir(parents=True, exist_ok=True)
            page.wait_for_selector(".requirement[data-status=undecided]")
            capture("mechanical")
        print("Dimensioned profile, native extrusion parameters, duplication and standard views passed.")

        page.locator("#studio-tab-pcb").click()
        assert page.locator("#mechanical-studio").is_hidden()
        assert page.locator("#pcb-canvas .pcb-pad").count() == 19
        page.locator("#pcb-route").click()
        page.locator('[data-pad="J1:1"]').click()
        page.locator('[data-pad="J1:2"]').click()
        assert "same assigned net" in page.locator("#pcb-status").inner_text()
        page.locator('[data-pad="U1:1"]').click()
        assert len(json.loads(page.evaluate("localStorage.getItem('emergence-pcb-v1')"))["design"]["routes"]) == 1
        page.locator("#pcb-minWidth").fill("1")
        page.locator("#pcb-minWidth").blur()
        assert page.locator('#pcb-checks [data-code="width"]').count() == 1
        page.locator("#pcb-undo").click()
        assert page.locator('#pcb-checks [data-code="width"]').count() == 0
        page.locator("#pcb-connectivity").click()
        assert "Logical net membership" in page.locator(".studio-netlist").inner_text()
        page.locator("#pcb-connectivity").click()
        page.locator("#pcb-route").click()  # Leave route mode before moving a component.
        page.locator("#pcb-studio .studio-tree-item").filter(has_text="J1 ·").click()
        page.locator("#pcb-x").fill("12")
        page.locator("#pcb-x").blur()
        assert len(json.loads(page.evaluate("localStorage.getItem('emergence-pcb-v1')"))["design"]["routes"]) == 0
        page.locator("#pcb-undo").click()
        with page.expect_download() as exported:
            page.locator("#pcb-export").click()
        board = json.loads(Path(exported.value.path()).read_text())
        assert len(board["design"]["routes"]) == 1
        corrupt = json.loads(json.dumps(board))
        corrupt["design"]["routes"][0]["to"] = "missing:1"
        page.locator("#pcb-import").set_input_files({"name": "bad.json", "mimeType": "application/json", "buffer": json.dumps(corrupt).encode()})
        page.wait_for_function("() => document.querySelector('#pcb-status').textContent.includes('Import rejected')")
        assert json.loads(page.evaluate("localStorage.getItem('emergence-pcb-v1')")) == board
        with page.expect_download() as bom:
            page.locator("#pcb-output").click()
        assert "Gate driver" in Path(bom.value.path()).read_text()
        page.reload(wait_until="networkidle")
        assert page.locator("#pcb-studio").is_visible()
        assert page.locator("#pcb-route-selection option").count() == 1
        # Valid import is undoable, and markup stays literal in both tree and SVG.
        imported = json.loads(json.dumps(board))
        imported["design"]["components"][0]["name"] = "<script>alert(1)</script>"
        page.locator("#pcb-import").set_input_files({"name": "board.json", "mimeType": "application/json", "buffer": json.dumps(imported).encode()})
        page.wait_for_function("() => document.querySelector('#pcb-status').textContent.includes('Imported design')")
        assert page.locator("#pcb-studio script").count() == 0
        assert "<script>" in page.locator("#pcb-studio .studio-tree-item").first.inner_text()
        page.locator("#pcb-undo").click()
        if screenshots:
            capture("pcb")
        print("PCB net selection, routing, width violations, move invalidation, undo, files, BOM and persistence passed.")

        page.locator("#studio-tab-silicon").click()
        page.locator("#silicon-studio .studio-tree-item").filter(has_text="SRAM macro").click()
        page.locator("#silicon-x").fill("400")
        page.locator("#silicon-x").blur()
        assert "Unlock" in page.locator("#silicon-status").inner_text()
        assert page.locator("#silicon-x").input_value() == "500"
        # Failed edits must leave the subsequent controls bound to the current document.
        page.locator("#silicon-fixed").uncheck()
        page.locator("#silicon-x").fill("200")
        page.locator("#silicon-x").blur()
        assert page.locator('#silicon-checks [data-code="overlap"]').count() == 1
        page.locator("#silicon-undo").click()
        assert page.locator('#silicon-checks [data-code="overlap"]').count() == 0
        page.locator("#silicon-studio .studio-tree-item").filter(has_text="CPU core").click()
        canvas = page.locator("#silicon-canvas")
        start = canvas.evaluate("""svg => {
            const p = svg.createSVGPoint(); p.x=200; p.y=220;
            const q=p.matrixTransform(svg.getScreenCTM()); return {x:q.x,y:q.y};
        }""")
        end = canvas.evaluate("""svg => {
            const p=svg.createSVGPoint(); p.x=240; p.y=240;
            const q=p.matrixTransform(svg.getScreenCTM()); return {x:q.x,y:q.y};
        }""")
        page.mouse.move(start["x"], start["y"])
        page.mouse.down()
        page.mouse.move(end["x"], end["y"], steps=5)
        page.mouse.up()
        assert page.locator("#silicon-x").input_value() == "140"
        assert page.locator("#silicon-y").input_value() == "120"
        page.locator("#silicon-undo").click()
        page.locator("#silicon-target").fill("10")
        page.locator("#silicon-target").blur()
        assert page.locator('#silicon-checks [data-code="density"]').count() == 1
        page.locator("#silicon-undo").click()
        with page.expect_download() as floorplan:
            page.locator("#silicon-output").click()
        svg = Path(floorplan.value.path()).read_text()
        assert "no PDK / no sign-off" in svg and "CPU core" in svg
        page.locator("#silicon-add").click()
        assert page.locator("#silicon-canvas [data-object]").count() == 5
        page.locator("#silicon-undo").click()
        if screenshots:
            capture("silicon")
        print("Floorplan locking, atomic failed edits, overlap, snap dragging, density, undo and SVG export passed.")

        page.locator("#studio-tab-silicon").focus()
        page.keyboard.press("Home")
        assert page.locator("#studio-tab-mechanical").get_attribute("aria-selected") == "true"
        page.keyboard.press("ArrowRight")
        assert page.locator("#studio-tab-pcb").get_attribute("aria-selected") == "true"
        for theme in ("light", "dark"):
            page.evaluate("(theme) => document.documentElement.dataset.theme=theme", theme)
            for width in (320, 390, 768, 1360):
                page.set_viewport_size({"width": width, "height": 844})
                for kind in ("mechanical", "pcb", "silicon"):
                    page.locator("#studio-tab-" + kind).click()
                    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth + 1"), (theme, width, kind)
        browser.close()
    assert not errors, errors
    print("Production CSP, keyboard workspace switching and responsive layouts in both themes passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url", nargs="?", default="http://127.0.0.1:18762/play/emergence/cad/")
    parser.add_argument("--screenshots", type=Path)
    args = parser.parse_args()
    check(args.url, args.screenshots)
