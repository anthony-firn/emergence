"""Mechanical must draw the project's populated PCB, with live edits and undo."""
import argparse
import hashlib
from pathlib import Path
from playwright.sync_api import sync_playwright
from test_cad_browser import POLICY


def check(url, screenshot=None):
    errors=[]
    with sync_playwright() as p:
        browser=p.chromium.launch()
        context=browser.new_context(viewport={'width':1440,'height':1000})
        def policy(route):
            response=route.fetch()
            route.fulfill(response=response,headers={**response.headers,'content-security-policy':POLICY})
        context.route('**/play/emergence/cad/',policy)
        page=context.new_page();page.on('pageerror',lambda error:errors.append(str(error)))
        page.goto(url+'?example=drone#mechanical',wait_until='networkidle')
        page.wait_for_function("() => window.emergence.solution && !document.body.classList.contains('solving')",timeout=60000)
        tree=page.locator('#feature-tree')
        assert 'Q1 · A high side' in tree.inner_text()
        assert 'Q6 · C low side' in tree.inner_text()
        assert 'Optical marker' not in tree.inner_text()
        assert 'thermal contact are not validated' in page.locator('#mechanical-assembly-note').inner_text()
        q6=page.evaluate("window.emergence.documents().pcb.components.find(c=>c.id==='Q6')")
        tree.get_by_role('button',name='▧ Q6 ·',exact=False).click()
        assert page.locator('#stl').is_disabled()
        page.locator('[data-camera=top]').click()
        canvas=page.locator('#assembly')
        before=hashlib.sha256(canvas.screenshot()).hexdigest()
        page.evaluate("async c => window.emergence.run('pcb.move',{id:c.id,at:[c.x+2,c.y]})",q6)
        after=hashlib.sha256(canvas.screenshot()).hexdigest()
        assert before!=after,'A PCB component move must redraw the Mechanical scene'
        page.locator('#command-undo').click()
        page.wait_for_function("c => window.emergence.documents().pcb.components.find(p=>p.id===c.id).x===c.x",arg=q6)
        assert hashlib.sha256(canvas.screenshot()).hexdigest()==before,'Undo must restore the drawn placement'
        page.locator('#exploded').check()
        assert hashlib.sha256(canvas.screenshot()).hexdigest()!=before
        page.locator('#exploded').uncheck()
        page.locator('[data-camera=iso]').click()
        if screenshot: canvas.screenshot(path=str(screenshot))
        for width in (320,390,1360):
            page.set_viewport_size({'width':width,'height':844})
            assert page.evaluate('document.documentElement.scrollWidth <= innerWidth+1'),width
        browser.close()
    assert not errors,errors
    print('Live PCB components, mechanical redraw, undo, exploded view, export guard, CSP and mobile widths passed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('url',nargs='?',default='http://127.0.0.1:8874/play/emergence/cad/')
    parser.add_argument('--screenshot',type=Path)
    args=parser.parse_args();check(args.url,args.screenshot)
