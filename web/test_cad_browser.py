"""Optional browser QA for the collaborative CAD demo using real browser peers."""
import argparse
import json
from pathlib import Path
from playwright.sync_api import sync_playwright

POLICY = ("default-src 'self'; script-src 'self' 'wasm-unsafe-eval'; "
          "style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; "
          "media-src 'self' data: blob:; worker-src 'self' blob:; connect-src 'self'; "
          "object-src 'none'; base-uri 'self'; frame-ancestors 'none'; form-action 'none'")


def check(url, screenshot):
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()

        def context():
            ctx = browser.new_context(viewport={"width": 1440, "height": 1100}, reduced_motion="reduce", color_scheme="dark", accept_downloads=True)

            def policy(route):
                response = route.fetch()
                route.fulfill(response=response, headers={**response.headers, "content-security-policy": POLICY})
            ctx.route("**/play/emergence/cad/", policy)
            ctx.on("page", lambda page: page.on("pageerror", lambda error: errors.append(str(error))))
            return ctx

        def open_page(ctx):
            page = ctx.new_page()
            page.goto(url, wait_until="networkidle")
            page.wait_for_selector(".requirement[data-status=undecided]")
            return page

        def domain(page, name):
            page.locator("#domains button").filter(has_text=name).click()

        def save(page, message):
            page.locator("#message").fill(message)
            page.locator("#commit").click()
            page.wait_for_function("() => document.querySelector('#commit').disabled")

        def revisions(page, count):
            page.wait_for_function("n => Number(document.querySelector('#revision-count').textContent) === n", arg=count)

        a_context, b_context = context(), context()
        a, b = open_page(a_context), open_page(b_context)
        a.locator("#author").fill("Ada mechanical")
        b.locator("#author").fill("Tess thermal")
        # Separate contexts do not share BroadcastChannel or localStorage.
        for page in (a, b):
            page.locator("#connect").click()
            page.locator("#peer-dialog summary").click()
            page.locator("#ice").fill("[]")  # Deterministic local candidates for this local peer test.
        a.locator("#make-offer").click()
        a.wait_for_function("() => document.querySelector('#peer-output').value.length > 100")
        b.locator("#peer-input").fill(a.locator("#peer-output").input_value())
        b.locator("#make-answer").click()
        b.wait_for_function("() => document.querySelector('#peer-output').value.length > 100")
        a.locator("#peer-input").fill(b.locator("#peer-output").input_value())
        a.locator("#finish-peer").click()
        for page in (a, b):
            page.wait_for_function("() => document.querySelector('#peer-count').textContent === '1'", timeout=30000)
            page.locator('[data-close="peer-dialog"]').click()
        print("Real WebRTC data channel connected between isolated contexts.")

        b.locator("#param-fins").fill("11")  # Unsaved draft on the reference revision.
        a.locator("#param-fins").fill("12")
        save(a, "Mechanical: twelve fins")
        revisions(b, 2)
        assert b.locator("#param-fins").input_value() == "11", "Incoming revision clobbered a draft"
        save(b, "Thermal: eleven fins")
        revisions(a, 3)
        a.wait_for_function("() => document.querySelector('#sync-status').textContent.includes('2 divergent')")
        a.locator("#tab-review").click()
        assert a.locator(".conflict").count() == 1
        assert a.locator("#merge").is_disabled()
        a.get_by_label("Resolve Heat-sink fins").select_option("source")
        a.wait_for_function("() => !document.querySelector('#merge').disabled")
        a.locator("#review-kind").select_option("approve")
        a.locator("#review-text").fill("Reviewed the pressure and thermal intervals. Keep eleven fins.")
        a.locator("#post-review").click()
        a.wait_for_selector(".note")
        b.wait_for_function("() => JSON.parse(localStorage.getItem('emergence-cad-project-v1')).notes.length === 1")
        # Exercise multiple RTC chunks whose JSON escaping expands packets.
        a.locator("#review-kind").select_option("comment")
        for _ in range(4):
            a.locator("#review-text").fill('"' * 2000)
            a.locator("#post-review").click()
            a.wait_for_function("() => document.querySelector('#review-text').value === ''")
        b.wait_for_function("() => JSON.parse(localStorage.getItem('emergence-cad-project-v1')).notes.length === 5")
        a.locator("#merge").click()
        revisions(a, 4)
        revisions(b, 4)
        b.wait_for_function("() => document.querySelector('#param-fins').value === '11'")
        print("Concurrent drafts, explicit conflict resolution, review sync and merge convergence passed.")

        a.locator("#tab-workspace").click()
        a.locator("#new-branch").click()
        a.locator("#branch-name").fill("thermal/cooling")
        a.locator('#branch-form button[type="submit"]').click()
        revisions(a, 5)
        a.locator("#param-velocity").fill("2")
        save(a, "Cooling: lower duct velocity")
        revisions(a, 6)
        a.locator("#branch").select_option(label="main")
        a.locator("#tab-review").click()
        a.wait_for_function("() => !document.querySelector('#merge').disabled")
        a.locator("#merge").click()
        revisions(a, 7)
        revisions(b, 7)
        a.locator("#tab-workspace").click()
        a.locator("#param-fins").fill("6")
        a.locator("#param-velocity").fill("3.5")
        a.wait_for_selector(".requirement[data-status=unavailable]")
        a.locator("#discard").click()
        a.wait_for_selector(".requirement[data-status=met]")
        domain(a, "Radiation")
        a.locator("#param-heat").fill("5")
        a.wait_for_function("() => document.querySelector('#model-status').textContent.includes('model runs')")
        assert "boundary input" in a.locator("#domain-description").inner_text()
        a.locator("#discard").click()
        domain(a, "Photonics")
        a.wait_for_function("() => document.querySelector('#optical-result').textContent.includes('reference discrepancy')")
        domain(a, "Mechanical")
        a.locator("#param-fins").fill("10")
        a.reload(wait_until="networkidle")
        assert a.locator("#param-fins").input_value() == "10", "Session draft was lost on reload"
        a.locator("#discard").click()
        print("Branch workflow, native model range failure, photonics and draft reload passed.")

        # Same-browser tabs use BroadcastChannel independently of WebRTC.
        b2 = open_page(b_context)
        b2.wait_for_function("() => Number(document.querySelector('#peer-count').textContent) >= 1")
        domain(b2, "Electrical")
        b2.locator("#param-current").fill("25")
        save(b2, "Electrical: derate current")
        revisions(b, 8)
        domain(b, "Electrical")
        b.wait_for_function("() => document.querySelector('#param-current').value === '25'")
        with b2.expect_download() as downloaded:
            b2.locator("#export").click()
        bundle = json.loads(Path(downloaded.value.path()).read_text())
        assert len(bundle["commits"]) == 8 and len(bundle["notes"]) == 5

        c = open_page(context())
        corrupt = json.loads(json.dumps(bundle))
        corrupt["commits"][-1]["params"]["current"] = 26
        c.locator("#import").set_input_files({"name": "tampered.json", "mimeType": "application/json", "buffer": json.dumps(corrupt).encode()})
        c.wait_for_function("() => document.querySelector('#toast').textContent.includes('verification')")
        assert c.locator("#revision-count").inner_text() == "1"
        c.locator("#import").set_input_files({"name": "project.json", "mimeType": "application/json", "buffer": json.dumps(bundle).encode()})
        revisions(c, 8)
        domain(c, "Electrical")
        assert c.locator("#param-current").input_value() == "25"
        domain(c, "Mechanical")
        with c.expect_download() as stl:
            c.locator("#stl").click()
        raw = Path(stl.value.path()).read_bytes()
        triangles = int.from_bytes(raw[80:84], "little")
        assert triangles > 100 and len(raw) == 84 + 50 * triangles
        assert b"millimetres" in raw[:80]
        assert c.locator(".support").get_attribute("href") == "/play/fairground/#/about"
        print(f"Local-tab sync, verified project round-trip and native STL export ({triangles} triangles) passed.")

        c.locator("#assembly").focus()
        c.keyboard.press("ArrowLeft")
        c.keyboard.press("+")
        c.locator("#tab-workspace").focus()
        c.keyboard.press("End")
        assert c.locator("#tab-history").get_attribute("aria-selected") == "true"
        c.keyboard.press("Home")
        assert c.locator("#tab-workspace").get_attribute("aria-selected") == "true"
        if screenshot:
            c.screenshot(path=str(screenshot), full_page=True)
        for theme in ("light", "dark"):
            c.emulate_media(color_scheme=theme)
            for width in (320, 390, 768):
                c.set_viewport_size({"width": width, "height": 844})
                for name in ("Mechanical", "Electrical", "Thermal", "Photonics", "Radiation"):
                    domain(c, name)
                    assert c.evaluate("document.documentElement.scrollWidth <= innerWidth + 1"), (theme, width, name)
                for tab in ("drawings", "review", "history", "workspace"):
                    c.locator(f"#tab-{tab}").click()
                    assert c.evaluate("document.documentElement.scrollWidth <= innerWidth + 1"), (theme, width, tab)
                c.locator("#connect").click()
                assert c.evaluate("document.querySelector('#peer-dialog').scrollWidth <= innerWidth"), (theme, width, "peers")
                c.locator('[data-close="peer-dialog"]').click()
        browser.close()
    assert not errors, errors
    print("Production CSP, keyboard controls and 320/390/768 px layouts in both themes passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url", nargs="?", default="http://127.0.0.1:18762/play/emergence/cad/")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    check(args.url, args.screenshot)
