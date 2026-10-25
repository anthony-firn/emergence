"""Optional end-to-end QA against the lab hosted inside evolvetech.dev.

Requires the Python Playwright package and its Chromium browser. No extra
dependencies are needed to build or use the lab itself.
"""
import argparse
from pathlib import Path
from playwright.sync_api import sync_playwright

POLICY = ("default-src 'self'; script-src 'self' 'wasm-unsafe-eval'; "
          "style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; "
          "media-src 'self' data: blob:; worker-src 'self' blob:; connect-src 'self'; "
          "object-src 'none'; base-uri 'self'; frame-ancestors 'none'; form-action 'none'")

def check(url, screenshot):
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=p.chromium.executable_path)
        page = browser.new_page(viewport={"width": 1400, "height": 1000}, reduced_motion="reduce")
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
        # Exercise the actual production policy even with a plain local server.
        def policy(route):
            response = route.fetch()
            route.fulfill(response=response, headers={**response.headers, "content-security-policy": POLICY})
        page.route("**/play/emergence/", policy)
        page.goto(url, wait_until="networkidle")
        page.wait_for_function("() => document.querySelector('#optics-result').textContent.includes('Rayleigh expansion')")
        page.select_option("#tolerance", "0.0000000001")
        page.wait_for_function("() => document.querySelector('#optics-result').textContent.includes('Unmet')")
        page.select_option("#tolerance", "0.00001")
        page.wait_for_function("() => !document.querySelector('#optics-result').textContent.includes('Unmet')")
        page.locator("#size").fill("1")
        page.wait_for_function("() => document.querySelector('#size-value').value === '10.00'")
        page.wait_for_function("() => document.querySelector('#optics-result').textContent.includes('Folded closure table')")

        page.click("#tab-flow")
        page.select_option("#flow-preset", "1")
        page.select_option("#resolution", "16")
        page.select_option("#reynolds", "40")
        page.wait_for_function("() => !document.querySelector('#play').disabled")
        assert float(page.locator("#flow-canvas").get_attribute("data-time")) == 0
        page.click("#play")
        page.wait_for_function("() => Number(document.querySelector('#flow-canvas').dataset.time) >= 0.12")
        page.click("#play")
        error = float(page.locator("#flow-canvas").get_attribute("data-error"))
        assert error < 1e-12, error
        print(f"Browser planar energy relative error: {error:.3e}")
        page.click("#reset")
        page.wait_for_function("() => document.querySelector('#flow-canvas').dataset.time === '0'")

        page.click("#tab-design")
        page.wait_for_selector(".check[data-status=undecided]")
        assert page.locator(".check[data-status=met]").count() == 2
        page.locator("#fins").fill("6")
        page.locator("#velocity").fill("3.5")
        page.click("#design-run")
        page.wait_for_selector(".check[data-status=unavailable]")
        assert "laminar" in page.locator("#design-results").inner_text().lower() or "Reynolds" in page.locator("#design-results").inner_text()

        # Keyboard tab navigation, both directions, including first/last.
        page.locator("#tab-design").focus()
        page.keyboard.press("Home")
        assert page.locator("#tab-optics").get_attribute("aria-selected") == "true"
        page.keyboard.press("End")
        assert page.locator("#tab-design").get_attribute("aria-selected") == "true"
        page.keyboard.press("ArrowRight")
        assert page.locator("#tab-optics").get_attribute("aria-selected") == "true"
        page.locator("#size").fill("-1")
        page.wait_for_function("() => document.querySelector('#optics-result').textContent.includes('Rayleigh expansion')")
        if screenshot:
            page.screenshot(path=str(screenshot))
        for theme in ("light", "dark"):
            page.emulate_media(color_scheme=theme)
            page.set_viewport_size({"width": 390, "height": 844})
            for tab in ("optics", "flow", "design"):
                page.click(f"#tab-{tab}")
                assert page.evaluate("document.documentElement.scrollWidth <= innerWidth + 1"), (theme, tab)
        page.goto(url + "#flow")
        page.wait_for_function("() => document.querySelector('#tab-flow').getAttribute('aria-selected') === 'true'")
        page.wait_for_function("() => !document.querySelector('#play').disabled")
        browser.close()
    assert not errors, errors
    print("Optics, flow, design, keyboard navigation, production CSP and mobile layouts passed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url", nargs="?", default="http://127.0.0.1:18761/play/emergence/")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    check(args.url, args.screenshot)
