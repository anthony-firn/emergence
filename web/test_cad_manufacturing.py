"""Authoring and manufacturing-output QA. Optional ezdxf independently audits DXF."""
import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
from playwright.sync_api import sync_playwright
from test_cad_browser import POLICY


def check(url):
    errors=[]
    with sync_playwright() as p:
        browser=p.chromium.launch()
        page=browser.new_page(viewport={"width":1500,"height":1100},color_scheme="dark",accept_downloads=True)
        page.on("pageerror",lambda e:errors.append(str(e)))
        def policy(route):
            response=route.fetch()
            route.fulfill(response=response,headers={**response.headers,"content-security-policy":POLICY})
        page.route("**/play/emergence/cad/",policy)
        page.goto(url+"#mechanical",wait_until="networkidle")
        page.wait_for_selector("#requirements .requirement[data-status=met]")
        page.locator("#param-width").fill("50")
        page.wait_for_function("() => document.querySelector('#handoff-values').textContent.includes('4.90 mm')")
        assert "50 × 50" in page.locator("#geometry-label").inner_text()
        page.locator("#add-part").click()
        page.locator("#part-name").fill("Insulating spacer")
        page.locator("#part-shape").select_option("cylinder")
        page.locator('#part-form button[type="submit"]').click()
        page.locator("#part-width").fill("14")
        page.locator("#part-width").blur()
        page.locator("#undo").click()
        assert page.locator("#part-width").input_value()=="12"
        page.locator("#redo").click()
        assert page.locator("#part-width").input_value()=="14"
        page.locator("#part-z").fill("5")
        page.locator("#part-z").blur()
        page.locator("#message").fill("Add a fourteen millimetre spacer and widen sink")
        page.locator("#commit").click()
        page.wait_for_function("() => document.querySelector('#commit').disabled")
        page.locator("#tab-drawings").click()
        assert page.locator("#drawing-part option:checked").inner_text()=="Insulating spacer"
        assert "⌀ 14" in page.locator("#drawing-preview").inner_text()
        for key,value in {"number":"ET-SP-001","material":"PEEK / unfilled","finish":"Ra 3.2 µm","tolerance":"0.1","prepared":"Designer","checked":"Shop checker"}.items():
            page.locator(f"#drawing-{key}").fill(value)
            page.locator(f"#drawing-{key}").blur()
        page.locator("#drawing-status").select_option("Released")
        assert "Save a revision" in page.locator("#drawing-readiness").inner_text()
        page.locator("#drawing-save").click()
        page.wait_for_function("() => document.querySelector('#drawing-readiness').textContent.includes('Declared release fields complete')")
        with page.expect_download() as dl:page.locator("#drawing-svg").click()
        svg=Path(dl.value.path()).read_text()
        root=ET.fromstring(svg)
        assert root.attrib["width"]=="420mm" and root.attrib["height"]=="297mm"
        assert "RELEASED BY NAMED CHECKER" in svg and "PEEK / unfilled" in svg and "GENERAL LINEAR TOLERANCE: ±0.1 mm" in svg
        assert "UNCOMMITTED DRAFT" not in svg
        with page.expect_download() as dl:page.locator("#drawing-dxf").click()
        dxf=Path(dl.value.path()).read_text()
        target=Path("/tmp/emergence-spacer-qa.dxf");target.write_text(dxf)
        assert "CIRCLE" in dxf and "40\n7\n" in dxf
        try:
            import ezdxf
        except ImportError:
            print("Optional independent DXF audit skipped: install ezdxf to enable it.")
        else:
            doc=ezdxf.readfile(target);audit=doc.audit()
            assert not audit.errors and not audit.fixes,([e.message for e in audit.errors],[f.message for f in audit.fixes])
            assert doc.units==4 and doc.modelspace().query("CIRCLE")[0].dxf.radius==7
        with page.expect_download() as dl:page.locator("#drawing-bom").click()
        rows=list(csv.reader(io.StringIO(Path(dl.value.path()).read_text())))
        assert len(rows)==3 and rows[-1][2]=="Insulating spacer" and rows[-1][3]=="PEEK / unfilled"
        page.pdf(path="/tmp/emergence-manufacturing-qa.pdf",prefer_css_page_size=True,print_background=True)
        pdf=Path("/tmp/emergence-manufacturing-qa.pdf").read_bytes();assert pdf.startswith(b"%PDF")
        info=subprocess.check_output(["pdfinfo","/tmp/emergence-manufacturing-qa.pdf"],text=True)
        assert "Pages:           1" in info and "(A3)" in info
        # Content edits cannot silently retain the old declared release.
        page.locator("#drawing-finish").fill("Machined, no coating")
        page.locator("#drawing-finish").blur()
        assert page.locator("#drawing-status").input_value()=="For review"
        page.locator("#tab-workspace").click()
        page.locator("#part-height").fill("8")
        page.locator("#part-height").blur()
        page.locator("#tab-drawings").click()
        assert "UNCOMMITTED DRAFT" in page.locator("#drawing-preview").inner_text()
        page.locator("#drawing-projection").select_option("first")
        assert "FIRST ANGLE PROJECTION" in page.locator("#drawing-preview").inner_text()
        page.locator("#drawing-notes").fill("<script>alert(1)</script> & keep literal markup")
        page.locator("#drawing-notes").blur()
        assert page.locator("#drawing-preview script").count()==0
        assert "<script>alert(1)</script>" in page.locator("#drawing-preview").inner_text()
        page.locator("#tab-workspace").click()
        with page.expect_download() as dl:page.locator("#stl").click()
        raw=Path(dl.value.path()).read_bytes();triangles=int.from_bytes(raw[80:84],"little")
        assert len(raw)==84+50*triangles and triangles>100
        for width in [320,390,768]:
            page.set_viewport_size({"width":width,"height":844})
            for tab in ["workspace","drawings","review","history"]:
                page.locator("#studio-tab-"+("history" if tab in ("review","history") else "mechanical")).click()
                page.locator(f"#tab-{tab}").click()
                assert page.evaluate("document.documentElement.scrollWidth<=innerWidth+1"),(width,tab)
        browser.close()
    assert not errors,errors
    print("Parametric authoring, undo/redo, revisioned shop drawings, release review, SVG/DXF/CSV/STL, A3 PDF and mobile layout passed.")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url",nargs="?",default="http://127.0.0.1:8080/play/emergence/cad/")
    check(parser.parse_args().url)
