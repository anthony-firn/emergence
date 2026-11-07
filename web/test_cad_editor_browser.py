"""Exercise the viewport editor under the production CSP, including real edits."""
import argparse
import hashlib
from pathlib import Path
from playwright.sync_api import sync_playwright
from test_cad_browser import POLICY


def check(url, screenshots=None):
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch(args=['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'])
        context = browser.new_context(viewport={'width': 1440, 'height': 1000}, accept_downloads=True)
        def policy(route):
            response = route.fetch()
            route.fulfill(response=response, headers={**response.headers, 'content-security-policy': POLICY})
        context.route('**/play/emergence/cad/', policy)
        page = context.new_page()
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.goto(url, wait_until='networkidle')
        page.wait_for_function("() => document.querySelector('#render-status').textContent.includes('triangles')", timeout=60000)
        assert page.locator('#studio-tab-render').get_attribute('aria-selected') == 'true'
        assert page.locator('#editor-tree [role=treeitem]').count() > 12
        assert page.evaluate('document.documentElement.scrollHeight === innerHeight')
        # Canvas and dock have real independent areas, rather than an overlay.
        canvas = page.locator('#render-canvas')
        assert canvas.bounding_box()['width'] > 950 and canvas.bounding_box()['height'] > 500
        assert page.locator('#editor-dock').bounding_box()['x'] >= canvas.bounding_box()['x'] + canvas.bounding_box()['width']

        # Select a real part, edit through the shared command registry, then undo.
        page.locator('#editor-tree').get_by_role('treeitem', name='⬡  Arm tubes', exact=True).click()
        assert page.locator('#editor-object-name').input_value() == 'Arm tubes'
        before = page.evaluate('window.emergence.documents().parts.parts[0]')
        page.locator('#editor-object-x').fill(str(before['x'] + 5))
        page.locator('#editor-object-x').blur()
        page.wait_for_function('(x) => window.emergence.documents().parts.parts[0].x === x', arg=before['x'] + 5)
        page.locator('#command-undo').click()
        page.wait_for_function('(x) => window.emergence.documents().parts.parts[0].x === x', arg=before['x'])
        # Collection placement must retain the original rotation, copies and phase.
        placement = page.evaluate('window.emergence.documents().project.scene.placements.parts')
        page.locator('#editor-x').fill(str(placement['at'][0] + 10))
        page.locator('#editor-x').blur()
        page.wait_for_function('(x) => window.emergence.documents().project.scene.placements.parts.at[0] === x', arg=placement['at'][0] + 10)
        after = page.evaluate('window.emergence.documents().project.scene.placements.parts')
        assert {k: v for k, v in placement.items() if k != 'at'} == {k: v for k, v in after.items() if k != 'at'}
        page.locator('#command-undo').click()
        page.wait_for_function('(p) => JSON.stringify(window.emergence.documents().project.scene.placements.parts) === JSON.stringify(p)', arg=placement)
        print('Part edits, scene placement, preserved rotations and undo passed.')

        page.locator('#editor-filter').fill('Arm tubes')
        assert page.locator('#editor-tree .editor-tree-children [role=treeitem]').count() == 1
        page.locator('#editor-filter').fill('no-such-object')
        assert page.locator('#editor-tree .editor-collection').count() == 0
        page.locator('#editor-filter').fill('')
        first = page.locator('#editor-tree [role=treeitem]').first
        first.focus(); page.keyboard.press('ArrowDown')
        assert page.evaluate("document.activeElement.closest('#editor-tree') !== null && !document.activeElement.classList.contains('editor-tree-root')")
        grip = page.get_by_role('separator', name='Resize properties')
        grip.focus(); page.keyboard.press('ArrowLeft')
        assert grip.get_attribute('aria-valuenow') == '320'
        page.locator('#editor-console-toggle').click()
        assert not page.locator('#command-input').is_visible()
        page.reload(wait_until='networkidle')
        assert grip.get_attribute('aria-valuenow') == '320'
        assert not page.locator('#command-input').is_visible()
        page.keyboard.press('Control+k')
        assert page.locator('#command-input').is_visible()
        assert page.locator('#command-input').evaluate('e => e === document.activeElement')
        page.locator('#editor-maximize').click()
        assert page.locator('#editor-area').evaluate("e => e.classList.contains('area-maximized')")
        page.keyboard.press('Escape')
        assert not page.locator('#editor-area').evaluate("e => e.classList.contains('area-maximized')")
        page.locator('#editor-select').focus(); page.keyboard.press('n')
        assert not page.locator('#editor-dock').is_visible()
        page.keyboard.press('n'); assert page.locator('#editor-dock').is_visible()
        # Typing shortcut letters must not toggle panels.
        page.locator('#editor-filter').fill('ng'); assert page.locator('#editor-dock').is_visible()
        page.locator('#editor-filter').fill('')
        page.keyboard.press('F3')
        assert page.locator('#commands-dialog').is_visible()
        assert page.locator('#commands-filter').evaluate('e => e === document.activeElement')
        page.keyboard.press('Escape')
        page.locator('#render-canvas').focus();page.keyboard.press('Control+Space')
        assert page.locator('#editor-area').evaluate("e => e.classList.contains('area-maximized')")
        page.keyboard.press('Escape')
        print('Filtering, keyboard navigation, pane resizing, saved layout and shortcuts passed.')

        page.wait_for_function("() => document.querySelector('#render-status').textContent.includes('triangles')", timeout=60000)
        views = []
        for name in ('top', 'front', 'right', 'iso'):
            page.locator('#editor-camera-' + name).click()
            views.append(hashlib.sha256(canvas.screenshot()).hexdigest())
        assert len(set(views)) == 4, 'Camera controls must produce distinct rendered views'
        # Pick through WebGL and show the picked collection in Properties.
        page.locator('#editor-camera-top').click()
        box = canvas.bounding_box()
        canvas.click(position={'x': box['width'] / 2, 'y': box['height'] / 2})
        assert page.locator('#editor-properties .panel-label').inner_text() != 'SCENE'
        page.locator('#editor-camera-iso').click()
        if screenshots:
            screenshots.mkdir(parents=True, exist_ok=True)
            page.screenshot(path=str(screenshots / 'emergence-editor-desktop.png'))
        # The compact layout exposes the spot plots immediately. Non-scaling
        # outlines need screen-sized dashes to keep the GPU responsive at focus.
        page.locator('#studio-tab-optics').click()
        z = page.evaluate('window.emergence.documents().optics.target.z')
        page.locator('#optics-focus').click()
        page.wait_for_function('(z) => window.emergence.documents().optics.target.z !== z', arg=z)
        page.locator('.spot').first.screenshot()
        page.locator('#command-undo').click()
        page.locator('#studio-tab-parts').click()
        page.wait_for_function("() => document.querySelector('#parts-mesh-info').dataset.mesh === 'ready'", timeout=60000)
        assert 'watertight' in page.locator('#parts-mesh-info').inner_text()
        print('Focused optics plots, undo, studio switching and part preview passed.')
        for theme in ('light', 'dark'):
            page.evaluate('(t) => document.documentElement.dataset.theme=t', theme)
            for width in (320, 390, 768, 1360):
                page.set_viewport_size({'width': width, 'height': 844})
                for studio in ('render', 'system', 'parts', 'mechanical', 'pcb'):
                    page.locator('#studio-tab-' + studio).click()
                    assert page.evaluate('document.documentElement.scrollWidth <= innerWidth + 1'), (theme, width, studio)
        if screenshots:
            page.set_viewport_size({'width': 390, 'height': 844})
            page.locator('#studio-tab-render').click()
            page.screenshot(path=str(screenshots / 'emergence-editor-mobile.png'), full_page=True)
        browser.close()
    assert not errors, errors
    print('Rendered camera views, WebGL selection, production CSP and responsive layouts passed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('url', nargs='?', default='http://127.0.0.1:8873/play/emergence/cad/')
    parser.add_argument('--screenshots', type=Path)
    args = parser.parse_args()
    check(args.url, args.screenshots)
