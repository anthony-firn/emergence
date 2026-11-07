// Desktop editor chrome. It is another client of the existing project bus:
// layout preferences stay local, while design edits keep the shared undo/history.
import { project, onSolved } from './app.js';
import { current } from './projects.js';
import { DOMAIN_INFO } from './domains/index.js';
import { MATERIALS } from './domains/parts.js';
import { node, button, numberField, textField, selectField } from './ui/kit.js';

export function mountEditor({ selectStudio, getStudio }) {
  const $ = id => document.getElementById(id);
  const prefsKey = 'emergence.editor.layout.v1';
  let prefs = {};
  try { prefs = JSON.parse(localStorage.getItem(prefsKey) || '{}') || {}; } catch {}
  if (typeof prefs !== 'object' || Array.isArray(prefs)) prefs = {};
  const save = () => { try { localStorage.setItem(prefsKey, JSON.stringify(prefs)); } catch {} };
  let selected = null, solved = null, pendingProperties = false;
  const header = document.querySelector('.app-header'), main = $('main');
  const scopeDialog = node('dialog'); scopeDialog.id = 'editor-about';
  const aboutTitle = node('div', undefined, 'dialog-heading');
  const heading = node('h2', 'Capabilities & design direction'); heading.id = 'editor-about-title';
  aboutTitle.append(heading, button('×', () => scopeDialog.close()));
  scopeDialog.setAttribute('aria-labelledby', heading.id);
  scopeDialog.append(aboutTitle, $('scope')); document.body.append(scopeDialog);
  document.querySelector('.app-footnote').remove();
  const menuBar = node('nav', undefined, 'editor-menus'); menuBar.setAttribute('aria-label', 'Editor menus');
  const menus = [];
  const command = id => () => { const target = $(id); if (target && !target.disabled) target.click(); };
  function menu(title, items) {
    const d = node('details', undefined, 'editor-menu'), summary = node('summary', title), list = node('div', undefined, 'editor-menu-items');
    for (const [name, fn, hint] of items) {
      const b = button(name, () => { d.open = false; fn(); }); if (hint) b.append(node('kbd', hint)); list.append(b);
    }
    d.append(summary, list); d.addEventListener('toggle', () => { if (d.open) for (const other of menus) if (other !== d) other.open = false; });
    menus.push(d); menuBar.append(d);
  }
  menu('File', [['Projects…', command('projects-open')], ['Rename project…', command('project-edit')], ['Open project file…', command('project-file')], ['Download project', command('project-download')], ['Export revision history', command('export')], ['Save revision…', command('global-save')]]);
  menu('Edit', [['Undo', () => project.bus.undo(), 'Ctrl Z'], ['Redo', () => project.bus.redo(), 'Ctrl Shift Z'], ['Command search…', () => { $('command-reference').click(); $('commands-filter').focus(); }, 'F3']]);
  const area = node('div', undefined, 'editor-area'); area.id = 'editor-area';
  const tools = node('aside', undefined, 'editor-tools'); tools.setAttribute('aria-label', 'Viewport tools');
  const centre = node('div', undefined, 'editor-workarea'); centre.id = 'editor-workarea';
  const dock = node('aside', undefined, 'editor-dock'); dock.id = 'editor-dock'; dock.setAttribute('aria-label', 'Scene outliner and properties');
  const grip = node('div', undefined, 'editor-resize'); grip.tabIndex = 0; grip.setAttribute('role', 'separator'); grip.setAttribute('aria-label', 'Resize properties'); grip.setAttribute('aria-orientation', 'vertical'); grip.setAttribute('aria-valuemin', '240'); grip.setAttribute('aria-valuemax', '480');
  const studioNodes = [...main.children].filter(e => e.matches('section[role=tabpanel]'));
  studioNodes[0].before(area); for (const n of studioNodes) centre.append(n);
  area.append(tools, centre, grip, dock);
  const projectbar = node('div', undefined, 'editor-projectbar'); document.querySelector('.studio-bar').before(projectbar);
  projectbar.append(document.querySelector('.repo-bar'), document.querySelector('.branch-bar'));
  const outliner = node('section', undefined, 'editor-outliner');
  const outHead = node('div', undefined, 'editor-pane-title'); outHead.append(node('strong', 'Outliner'), node('span', 'Scene collection', 'editor-pane-kind'));
  const filter = node('input'); filter.type = 'search'; filter.id = 'editor-filter'; filter.placeholder = 'Search collection'; filter.setAttribute('aria-label', 'Search scene collection');
  const tree = node('div', undefined, 'editor-tree'); tree.id = 'editor-tree'; tree.setAttribute('role', 'tree'); tree.setAttribute('aria-label', 'Project collection');
  outliner.append(outHead, filter, tree);
  const propHead = node('div', undefined, 'editor-pane-title'); propHead.append(node('strong', 'Properties'), node('span', 'Context', 'editor-pane-kind'));
  const inspector = node('section', undefined, 'editor-properties'); inspector.id = 'editor-properties'; inspector.setAttribute('aria-label', 'Selection properties');
  dock.append(outliner, propHead, inspector);

  const consoleDock = $('command-dock'); consoleDock.classList.add('editor-console');
  const consoleHead = node('div', undefined, 'editor-console-heading');
  const toggleConsole = button('⌃ Console', () => setConsole(!consoleDock.classList.contains('is-collapsed'))); toggleConsole.id = 'editor-console-toggle';
  const history = button('Revisions & reviews', () => selectStudio('history'));
  const activity = node('span', 'No commands in this session', 'editor-activity'); activity.id = 'editor-activity';
  consoleHead.append(toggleConsole, history, activity); consoleDock.prepend(consoleHead);
  const statusbar = node('footer', undefined, 'editor-statusbar'); statusbar.id = 'editor-statusbar';
  const navigationHint = node('span', 'Drag orbit · Shift drag pan · Wheel zoom · F3 commands');
  const solverStatus = node('span', 'Rust / WASM · solving…'); solverStatus.id = 'editor-solver-state';
  statusbar.append(navigationHint, $('sync-status'), solverStatus); consoleDock.after(statusbar);
  function setConsole(collapsed) {
    consoleDock.classList.toggle('is-collapsed', collapsed); toggleConsole.textContent = collapsed ? '⌄ Console' : '⌃ Console';
    toggleConsole.setAttribute('aria-expanded', String(!collapsed));
    prefs.consoleCollapsed = collapsed; save();
  }
  function setProperties(hidden) { area.classList.toggle('properties-hidden', hidden); prefs.propertiesHidden = hidden; save(); }
  function maximize() { area.classList.toggle('area-maximized'); document.dispatchEvent(new Event('emergence:editor-resize')); }
  function resetLayout() { prefs = {}; area.classList.remove('properties-hidden', 'area-maximized'); setWidth(300); setConsole(false); filter.value = ''; renderTree(); }
  menu('View', [['Maximize editor', maximize, 'Ctrl Space'], ['Toggle properties', () => setProperties(!area.classList.contains('properties-hidden')), 'N'], ['Toggle console', () => setConsole(!consoleDock.classList.contains('is-collapsed'))], ['Reset layout', resetLayout]]);
  menu('Help', [['Project tour', command('start-tour')], ['Command reference', command('command-reference')], ['Capabilities & limits…', () => scopeDialog.showModal()], ['Physics lab', () => location.assign('../')]]);
  header.querySelector('.app-brand').after(menuBar);
  document.addEventListener('pointerdown', e => { if (!menuBar.contains(e.target)) for (const menu of menus) menu.open = false; });
  document.addEventListener('keydown', e => { if (e.key === 'Escape') { for (const menu of menus) menu.open = false; if (area.classList.contains('area-maximized')) area.classList.remove('area-maximized'); } });

  const viewport = (action, view) => {
    if (getStudio() !== 'render') selectStudio('render');
    document.dispatchEvent(new CustomEvent('emergence:viewport', { detail: { action, view } }));
  };
  function tool(glyph, label, fn, id) { const b = button(glyph, fn); b.title = label; b.setAttribute('aria-label', label); if (id) b.id = id; tools.append(b); return b; }
  tool('↖', 'Inspect the 3D scene', () => { if (getStudio() !== 'render') selectStudio('render'); $('render-canvas')?.focus(); }, 'editor-select');
  tool('✥', 'Edit collection placement in Properties', () => { setProperties(false); if (!['parts', 'mech', 'pcb', 'silicon', 'optics', 'flight'].includes(selected?.domain)) selected = { domain: project.active().find(d => ['parts', 'mech', 'pcb', 'silicon', 'optics', 'flight'].includes(d)) }; renderProperties(); $('editor-x')?.focus(); }, 'editor-transform');
  tools.append(node('span', '', 'editor-tool-divider'));
  tool('⌖', 'Fit scene', () => viewport('fit'), 'editor-fit');
  for (const [glyph, name] of [['Z', 'top'], ['Y', 'front'], ['X', 'right'], ['◇', 'iso']]) tool(glyph, `${name[0].toUpperCase() + name.slice(1)} view`, () => viewport('camera', name), `editor-camera-${name}`);
  tools.append(node('span', '', 'editor-tool-divider'));
  tool('◌', 'Toggle X-ray in Layout', () => { if (getStudio() !== 'render') selectStudio('render'); $('render-xray')?.click(); }, 'editor-xray');
  tool('⤢', 'Maximize editor', maximize, 'editor-maximize');
  function setWidth(value) { const width = Math.max(240, Math.min(480, Number(value) || 300)); area.style.setProperty('--editor-dock-width', width + 'px'); grip.setAttribute('aria-valuenow', String(width)); prefs.dockWidth = width; save(); }
  let resize = null;
  grip.addEventListener('pointerdown', e => { resize = e.clientX; grip.setPointerCapture(e.pointerId); e.preventDefault(); });
  grip.addEventListener('pointermove', e => { if (resize !== null) { setWidth((prefs.dockWidth || 300) + resize - e.clientX); resize = e.clientX; } });
  grip.addEventListener('pointerup', () => resize = null); grip.addEventListener('pointercancel', () => resize = null);
  grip.addEventListener('keydown', e => { if (['ArrowLeft', 'ArrowRight', 'Home'].includes(e.key)) { e.preventDefault(); setWidth(e.key === 'Home' ? 300 : (prefs.dockWidth || 300) + (e.key === 'ArrowLeft' ? 20 : -20)); } });
  function choose(selection) {
    const inTree = tree.contains(document.activeElement);
    selected = selection; renderTree(); renderProperties(true);
    if (inTree) tree.querySelector('[aria-selected=true]')?.focus();
    document.dispatchEvent(new CustomEvent('emergence:selection-view', { detail: selection }));
  }
  function items(domain) {
    const docs = project.store.docs, d = docs[domain];
    if (domain === 'parts') return (d?.parts || []).map(p => ({ id: p.id, name: p.name, icon: '⬡' }));
    if (domain === 'pcb') return (d?.components || []).map(c => ({ id: c.id, name: `${c.id} · ${c.name}`, icon: '▦' }));
    if (domain === 'silicon') return (d?.blocks || []).map(c => ({ id: c.id, name: c.name, icon: '◈' }));
    if (domain === 'architecture') return (d?.rooms || []).map(c => ({ id: c.id, name: c.name, icon: '⌂' }));
    return [];
  }
  function renderTree() {
    const query = filter.value.trim().toLowerCase(); tree.replaceChildren();
    const root = button('▾  ' + project.store.docs.project.name, () => choose(null), 'editor-tree-root'); root.setAttribute('role', 'treeitem'); root.setAttribute('aria-selected', String(!selected)); tree.append(root);
    for (const domain of project.active()) {
      const info = DOMAIN_INFO[domain], children = items(domain), matches = !query || info.title.toLowerCase().includes(query);
      const visible = matches ? children : children.filter(c => c.name.toLowerCase().includes(query));
      if (!matches && !visible.length) continue;
      const group = node('details', undefined, 'editor-collection'); group.open = !!query || selected?.domain === domain || domain === 'parts';
      const summary = node('summary');
      const select = button(`${info.icon || '◇'}  ${info.title}`, () => choose({ domain }), 'editor-tree-item'); select.setAttribute('role', 'treeitem'); select.setAttribute('aria-selected', String(selected?.domain === domain && !selected?.id));
      const count = node('span', String(children.length || '·'), 'editor-item-count'); summary.append(select, count); group.append(summary);
      const nodes = node('div', undefined, 'editor-tree-children'); nodes.setAttribute('role', 'group');
      for (const item of visible) {
        const b = button(`${item.icon}  ${item.name}`, () => choose({ domain, id: item.id, name: item.name }), 'editor-tree-item'); b.setAttribute('role', 'treeitem'); b.setAttribute('aria-selected', String(selected?.domain === domain && selected?.id === item.id)); nodes.append(b);
      }
      if (!children.length) nodes.append(button('Open ' + info.title + ' workspace →', () => selectStudio(info.studio), 'editor-tree-open'));
      group.append(nodes); tree.append(group);
    }
  }
  function renderProperties(force = false) {
    if (!force && inspector.contains(document.activeElement) && document.activeElement.matches('input,select,textarea')) { pendingProperties = true; return; }
    pendingProperties = false;
    inspector.replaceChildren();
    const docs = project.store.docs, d = selected?.domain, info = DOMAIN_INFO[d];
    const status = node('p', '', 'studio-status'); status.id = 'editor-property-status';
    if (!d || !project.active().includes(d)) {
      inspector.append(node('p', 'SCENE', 'panel-label'), node('h2', docs.project.name), node('p', docs.project.summary || 'Select a collection or an object to inspect it.', 'small muted'));
      const facts = node('dl', undefined, 'editor-facts');
      for (const [label, value] of [['Domains', project.active().length], ['Units', 'mm · kg · s'], ['Storage', 'This device'], ['Project', current.meta.source?.example ? 'Editable example' : 'Local project']]) facts.append(node('dt', label), node('dd', String(value)));
      inspector.append(facts, button('Edit project details…', command('project-edit')), button('Add a domain…', command('studio-add')));
      return;
    }
    inspector.append(node('p', info.title.toUpperCase(), 'panel-label'), node('h2', selected.name || info.title));
    const part = d === 'parts' && docs.parts.parts.find(p => p.id === selected.id);
    if (part) {
      textField(inspector, { label: 'Name', id: 'editor-object-name', value: part.name, op: 'parts.update', args: name => ({ id: part.id, name }), status });
      selectField(inspector, { label: 'Material', id: 'editor-material', value: part.material, options: Object.entries(MATERIALS).map(([id, m]) => [id, m.title]), op: 'parts.update', args: material => ({ id: part.id, material }), status });
      inspector.append(node('h3', 'Object transform'));
      for (const axis of ['x', 'y', 'z']) numberField(inspector, { label: axis.toUpperCase(), id: 'editor-object-' + axis, value: part[axis], unit: 'mm', min: -2000, max: 2000, op: 'parts.update', args: value => ({ id: part.id, [axis]: value }), status });
    }
    const placement = docs.project.scene.placements[d];
    if (['parts', 'mech', 'pcb', 'silicon', 'optics', 'flight'].includes(d)) {
      inspector.append(node('h3', 'Collection placement'), node('p', `Positions the ${info.title} collection in the shared 3D scene.`, 'small muted'));
      const field = (axis, value) => {
        const label = node('label', axis.toUpperCase()), input = node('input'); input.type = 'number'; input.step = 'any'; input.value = value ?? ''; input.id = 'editor-' + axis; label.htmlFor = input.id; label.append(input); inspector.append(label);
        input.addEventListener('change', async () => {
          if (!input.checkValidity()) return;
          const p = project.store.docs.project.scene.placements[d], at = [...(p?.at || [0, 0, 0])]; at[['x', 'y', 'z'].indexOf(axis)] = input.value === '' ? 0 : Number(input.value);
          try { await project.bus.run('project.place', { domain: d, x: at[0], y: at[1], z: at[2], rx: p?.rx, rz: p?.rz, phase: p?.phase, ...(p?.count !== undefined ? { count: String(p.count) } : {}) }); }
          catch (error) { status.textContent = error.message; }
        });
      };
      for (const [i, axis] of ['x', 'y', 'z'].entries()) field(axis, placement?.at[i]);
      if (!placement) inspector.append(node('p', 'Automatic layout until a placement is entered.', 'small muted'));
    }
    inspector.append(button('Open ' + info.title + ' workspace →', () => {
      selectStudio(info.studio);
      if (part) document.querySelectorAll('#parts-studio .parts-list .studio-tree-item')[docs.parts.parts.findIndex(p => p.id === part.id)]?.click();
    }), status);
    const checks = solved?.checks.filter(c => c.domain === d) || [];
    if (checks.length) {
      inspector.append(node('h3', 'Domain checks'));
      const list = node('ul', undefined, 'editor-checks');
      for (const c of checks.slice(0, 8)) { const item = node('li', `${c.status === 'met' ? '✓' : c.status === 'unmet' ? '×' : '·'} ${c.text}`); item.dataset.status = c.status; list.append(item); } inspector.append(list);
    }
  }
  inspector.addEventListener('focusout', () => { setTimeout(() => { if (pendingProperties) renderProperties(); }, 0); });
  tree.addEventListener('keydown', e => {
    const rows = [...tree.querySelectorAll('[role=treeitem]')].filter(n => n.getClientRects().length);
    const i = rows.indexOf(document.activeElement);
    if (['ArrowDown', 'ArrowUp', 'Home', 'End'].includes(e.key)) {
      e.preventDefault(); rows[e.key === 'Home' ? 0 : e.key === 'End' ? rows.length - 1 : Math.max(0, Math.min(rows.length - 1, i + (e.key === 'ArrowDown' ? 1 : -1)))]?.focus();
    } else if (['ArrowLeft', 'ArrowRight'].includes(e.key)) {
      const group = document.activeElement.closest('.editor-collection');
      if (group) { e.preventDefault(); group.open = e.key === 'ArrowRight'; if (!group.open) group.querySelector('[role=treeitem]')?.focus(); }
    }
  });
  filter.addEventListener('input', renderTree);
  document.addEventListener('emergence:selection', e => { selected = e.detail; renderTree(); renderProperties(true); });
  document.addEventListener('emergence:studio', e => {
    document.body.dataset.editorStudio = e.detail;
    const tab = document.querySelector(`[data-studio="${e.detail}"]`); tab?.scrollIntoView({ block: 'nearest', inline: 'nearest' });
  });
  project.store.on(({ ids }) => { if (selected?.domain && !project.active().includes(selected.domain)) selected = null; renderTree(); renderProperties(); });
  onSolved(s => { solved = s; const failures = s.checks.filter(c => ['unmet', 'unavailable'].includes(c.status)).length; solverStatus.textContent = `${s.converged ? 'Solved' : 'Not settled'} · ${s.iterations} passes · ${failures} failing / unavailable checks`; renderProperties(); });
  project.bus.on(event => { activity.textContent = event.type === 'error' ? `${event.op}: ${event.error.message}` : `${event.type}: ${event.op || event.label || 'project updated'}`; });
  document.addEventListener('keydown', e => {
    const editing = e.target.closest('input,textarea,select,[contenteditable=true]');
    if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'k') { e.preventDefault(); setConsole(false); $('command-input').focus(); $('command-input').select(); return; }
    if (e.key === 'F3') { e.preventDefault(); setConsole(false); $('command-reference').click(); $('commands-filter').focus(); }
    if (editing || document.querySelector('dialog[open]')) return;
    if (e.ctrlKey && e.code === 'Space') { e.preventDefault(); maximize(); }
    else if (!e.ctrlKey && !e.metaKey && !e.altKey && e.key.toLowerCase() === 'n') { e.preventDefault(); setProperties(!area.classList.contains('properties-hidden')); }
    else if (!e.ctrlKey && !e.metaKey && !e.altKey && e.key.toLowerCase() === 'g') { e.preventDefault(); $('editor-transform').click(); }
  });
  setWidth(prefs.dockWidth || 300); setConsole(!!prefs.consoleCollapsed); setProperties(!!prefs.propertiesHidden);
  renderTree(); renderProperties();
  document.body.dataset.editorStudio = getStudio();
  if (!location.hash) selectStudio('render');
}
