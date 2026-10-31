# Collaborative CAD and manufacturing handoff

`/play/emergence/cad/` is a browser engineering workspace built on Emergence's
native CAD kernel and project graph. Users can edit the reference cooling
assembly, add parametric parts, branch, review, prepare drawings and exchange
revisions. It is an early mechanical authoring and handoff tool, not a qualified
replacement for all CAD/CAM and manufacturing workflows.

## Browser editors and professional workflow direction

The workspace is now one cross-domain project: an inspection quadcopter whose
flight, mechanical, parts, PCB, silicon, optics, fluids and architecture
documents are linked by unit-checked ports, solved together, and saved as
content-addressed blobs in every revision. See `analysis/cross_domain_cad.md`
for the architecture, the public reference designs and every check. The
sections below describe the mechanical workflow and the earlier PCB and silicon
editors, which are now studios of that project.

The website now exposes three workspaces, with direct links at
cad/#mechanical, cad/#pcb and cad/#silicon. Their interaction patterns
take cues from the [SolidWorks feature tree](https://help.solidworks.com/2023/english/Solidworks/sldworks/c_FeatureManager_Design_Tree_Conventions.htm),
[Altium interactive routing](https://www.altium.com/documentation/altium-designer/pcb/routing/interactive)
and [Cadence Innovus implementation](https://www.cadence.com/en_US/home/resources/datasheets/innovus-implementation-system-ds.html).
These references describe the direction, not equivalent capability.

Mechanical authoring adds an origin-centred rectangular/circular XY profile
with explicit dimensions and extrusion distance, a profile preview with exact
analytic volume, feature descriptions, duplicate extrusion, standard orthographic
views and rebuild/check commands. The resulting parts remain the same native
box/cylinder solids, drawing definitions and revisioned features. This does not
introduce a general geometric constraint solver, face references or B-rep.

The PCB editor has editable illustrative through-hole footprints, 90-degree
rotations, pad-net assignments, a logical connectivity list, top/bottom copper
visibility, airwires and explicit pad-to-pad connections with a 45-degree bend.
Geometry checks cover circular pad and straight track clearance, track width,
copper-to-edge distance and body containment. Crossings between different nets
on one layer are detected; crossings on different layers do not imply a short.
Moving, rotating or reassigning a component removes its attached connections.
Connectivity counts only centre-terminated explicit pad links; incidental
same-net copper intersections are not inferred. There is no obstacle avoidance,
schematic circuit simulation, via/plane model or Gerber/drill output. Component
sizes and pad radii are illustrative estimates, not qualified land patterns.

The silicon editor places rectangular functional blocks in micrometres, with
snap, editable size/position, locking, core keepout, rectangular spacing halos
and a requested block-area target. It detects overlapping blocks and core-boundary
violations. Area usage is the **sum of requested block areas**, divided by the
inner core area; overlaps are reported separately and are not subtracted from
that sum. SVG export is an annotated concept floorplan. No PDK, technology
mapping, cell placement/routing, extraction, timing, DRC/LVS or GDSII is provided
by this browser editor. The separate native chip flow retains its existing
synthesis, equivalence, timing and power validation in analysis/chip.md.

PCB and silicon keep their emergence-studio/1 JSON import/export, and are now
domains of the shared project: their documents are revisioned, merged by object
id and synchronised with peers along with every other domain. Failed imports and edits leave the document intact.
Unavailable browser storage is reported, with file export still available.

The studio model tests check exact segment-distance cases (crossing, parallel,
collinear, degenerate and endpoint distance), footprint rotation, connectivity,
same-layer copper crossings, pad clearance and rectangular area/bounds/halos.
The retained exact distances 0, 3 and 5 and the 12.5% area fixture have zero
absolute discrepancy. The circular clearance fixture accepts exactly 0.4 mm
and detects a 0.001 mm deficit. These are geometric algorithm checks, not a
physical accuracy certificate or fabrication validation.

Run node --test web/test_cad_studio.mjs. The browser suite
web/test_cad_studios_browser.py exercises profile dimensions, routing, move
invalidation, locking, snapped dragging, rejected edits, persistence,
JSON/BOM/SVG, production CSP, keyboard tabs and both themes at 320–1360 px.

## Geometry and cross-domain models

The reference sink has editable width, length, base thickness, fin count,
height and thickness. Those dimensions drive the viewport, analytic drawing,
DXF profile, native STL and the actual motor-driver thermal channel model.
Current, trace width, inlet temperature and duct velocity feed the original
unit-checked electrothermal graph. Model failures, intervals crossing limits
and unsupported flow regimes remain visible.

Prescribed optical power times absorbed fraction contributes deposited watts
at Q1, as does the external heat boundary. These watts are added to switching
loss, retaining its conservative 30% uncertainty interval. This is a
conservation-based boundary model, not a resolved absorption or nuclear
transport calculation. The separate Mie model at 550 nm does not determine
that absorption fraction. Constants and their sources remain in the original
native model code.

The feature tree supports up to sixteen additional box or cylinder parts,
with stable IDs, names, disciplines, dimensions and placement. Select parts in
the tree or viewport, show/hide them locally, edit, delete and undo/redo drafts.
Custom parts are geometric solids; they do not silently become new thermal,
optical or structural couplings. Power-board components are schematic context.
The presets are editable design starting points, not qualified devices.

STL export runs the real `cad::Solid` surface-net mesher. Sink cell spacing
resolves both fin thickness and channel gap. Exports requiring a cell below
0.2 mm or more than two million grid nodes report the browser limit instead
of claiming to preserve those details. Meshing rounds sharp edges. Dimensions
and positions are in millimetres; custom-part STL retains assembly placement.
A native test checks reference mesh volume against the exact sum of solids
and watertight topology; another checks a cylinder against pi r² h.

## Mechanical shop drawings

The Shop drawings view supports the finned extrusion, boxes and cylinders.
Its analytic dimensions come from the design snapshot, never the rounded mesh.
It produces:

- An A3 landscape vector sheet with aligned top/front/right views, visible and
  centre lines where appropriate, dimension and extension lines, arrowheads,
  repeated fin details, projection symbol and title block.
- First- or third-angle projection. Written dimensions are millimetres.
- Part-specific number, title, revision, material/grade/condition, finish,
  general linear tolerance, preparer, checker, notes and declared sheet status.
- SVG and browser print/PDF preserving the vector sheet and physical A3 size.
- A true-size AutoCAD 2000 DXF profile with millimetre insertion units: the
  extrusion end section for a sink; the local top profile for a box/cylinder.
- A mechanical parts CSV listing the supported solids, drawing identifiers,
  material fields, quantities and content revision. It excludes schematic
  board components for which no fabrication definition exists.

Drawings start as drafts. A declared Released sheet requires filled material,
finish, tolerance, drawing number, title, revision, preparer and checker fields.
An unsaved design is stamped UNCOMMITTED DRAFT. Geometry edits return released
sheets to For review; edits to released drawing fields do likewise. Undo may
restore the exact earlier snapshot, including its status. Sheet declarations
are entered by the team; the prototype does not verify the checker identity.

Drawing metadata is included in content-addressed commits and merge review.
A sheet carries the relevant content hash. Its hash and declared revision are
traceability aids, not proof that an engineer approved the design. General
linear tolerance is an explicit user-entered value; the software does not
invent process capability, glass specifications, fit classes or GD&T controls.

The projection and drawing conventions are informed by
[ISO 5456-2](https://www.iso.org/obp/ui#iso:std:iso:5456:-2:ed-1:v1:en).
[ASME Y14.5](https://www.asme.org/codes-standards/find-codes-standards/dimensioning-and-tolerancing)
and the [ASME Y14 series](https://www.asme.org/codes-standards/y14-standards)
cover additional definition and drawing practices. This prototype does not
claim conformance certification or full implementation of those standards.
GD&T, arbitrary sketches and constraints, STEP/B-rep editing, machining setups,
CAM toolpaths, sheet-metal unfolding, optical fabrication sheets and PCB
Gerber/drill packages remain unfinished. DXF is geometry for downstream CAD/CAM,
not a machine program. Independent CAD/CAM and shop review is still required.

## Revisions and compatibility

`web/cad/repository.js` defines `emergence-cad/1`. Each immutable revision
contains parent hashes, branch, parameter snapshot, author, discipline, message
and timestamp. New revisions also carry custom parts and drawing metadata.
SHA-256 covers canonical JSON with recursively sorted object keys.

The original six-parameter root and legacy content hashes remain unchanged.
Missing expanded parameters acquire the reference defaults for evaluation and
diffing, without rewriting old records. Exports from the first prototype can
therefore be imported. Both collaborating clients should run the current app;
old clients do not understand the additional fields.

Replicas union verified revisions and review notes. Branch tips derive from
ancestry, so concurrent edits produce multiple heads instead of overwriting
one another. Three-way merge uses the unique best common ancestor. Independent
parameter and part changes merge; conflicting parameter changes, edit/delete
part conflicts, and conflicting sheet metadata require an explicit choice.
Criss-cross histories with multiple best bases stop for reconciliation.
Merged electrothermal inputs are evaluated before enabling merge. Requirement
failures remain visible; recording a revision does not certify the device.
Merge commits preserve both parents and retain the target as first parent.

Incoming fast-forwards update clean workspaces. Unsaved drafts stay on their
original base. Committed history is in localStorage; active revision and draft
are in sessionStorage, with fifty in-memory undo states. Browser storage can
be cleared. Export includes committed history and reviews, not drafts. There
is no cloud backup. Names and domain ownership are self-reported; hashes check
content integrity, not authorship. No signed identities or enforced permissions
are implemented yet.

Imports validate schemas, IDs, hashes, supported shapes, bounded coordinates,
drawing fields and connected ancestry before changing the repository. The
combined compact project limit is 3 MB, 512 revisions and 1,000 notes.

## Peer connections

Same-origin browser tabs exchange revisions with BroadcastChannel. Separate
browsers use ordered WebRTC data channels with manual offer/answer pairing.
Cloudflare's public STUN service is the default network discovery aid; the
user can substitute ICE server JSON or use an empty array for local candidates.
Some networks need user-supplied TURN relaying. No TURN credentials, signaling
server, permanent discovery, offline mailbox or background replication is
provided. See [Cloudflare's STUN endpoint](https://developers.cloudflare.com/realtime/sfu/get-started/connection-patterns/)
and [WebRTC connectivity](https://developer.mozilla.org/en-US/docs/Web/API/WebRTC_API/Connectivity).

Peers exchange bounded full snapshots, chunked in 12,000-character segments
for RTC transport. Imports serialize validation. Full history travels
transitively across connected peers. Everyone invited receives the whole
project and review notes. Pair only intended collaborators. Keep a project
export if it must survive all peers closing their browsers.

## Funding

Support reuses `/play/fairground/#/about` and its existing donation flow. The
public explorer link follows Fairground’s current development address:
`WMQHGXHHWQVUN2VIZQJ5OO2GSH5Q33CVABJUCLJHP4UHLTTMS2HE3W6FEM`.
The recipient is defined by `src/defaults.ts` in Fairground. CAD does not hold wallet keys or implement a separate checkout.

## Validation

From the repository root, with the built site served locally:

```sh
(cd fluid && cargo test -p emergence-web)
node --test web/test_cad_repository.mjs
python3 web/test_cad_browser.py http://127.0.0.1:8080/play/emergence/cad/
python3 web/test_cad_manufacturing.py http://127.0.0.1:8080/play/emergence/cad/
```

Browser QA uses optional Python Playwright and Chromium; manufacturing QA also
uses `pdfinfo`. Optional `ezdxf` independently loads and audits the generated
DXF, checks insertion units and verifies the actual profile geometry. These
are development tools, not application dependencies.

Checks cover optical-to-thermal power equivalence, geometry propagation,
analytic section area, primitive volume/topology, legacy histories, concurrent
parts and sheet conflicts, real RTC peers in isolated browser contexts,
draft preservation, large escaped review transfers, export/import, tamper
rejection, authoring and undo/redo, declared release review, SVG escaping,
DXF/CSV/STL downloads, a single A3 vector PDF, production CSP, keyboard controls
and small-screen layouts. These checks establish the stated prototype
behaviour; they do not establish compatibility with every shop or CAM package.
