//! Anisotropic cells and block culling: the two things that make a thin
//! feature in a wide model meshable.

use std::f64::consts::PI;

use cad::kits::pcb::{Board, Feature, PcbLayer, OZ_COPPER};
use cad::meshing::{mesh_solid, mesh_solid_anisotropic, surface_nets, surface_nets_sparse, Grid3};
use cad::sketch::Region2;
use cad::solid::Solid;

#[test]
fn block_culling_gives_exactly_the_dense_mesh() {
    // The claim is identity, not approximation: a 1-Lipschitz field cannot
    // change sign inside a block whose centre is farther from zero than
    // the block's own half-diagonal, so culling such blocks cannot lose a
    // surface. Cells are visited in the dense order, so even the vertex
    // numbering matches.
    let solid = Solid::box_between([-0.2, -0.12, -0.02], [0.2, 0.12, 0.02])
        .minus(Solid::cylinder([-0.1, 0.0, -0.1], [-0.1, 0.0, 0.1], 0.03))
        .minus(Solid::cylinder([0.1, 0.0, -0.1], [0.1, 0.0, 0.1], 0.03))
        .plus(Solid::sphere([0.15, 0.08, 0.02], 0.03));
    let grid = Grid3::covering(solid.bbox(), 0.004, 2).unwrap();

    let dense = surface_nets(&solid, &grid, 2);
    for block in [2usize, 4, 8, 16] {
        let (sparse, stats) = surface_nets_sparse(&solid, &grid, 2, block);
        assert_eq!(sparse.vertices.len(), dense.vertices.len(), "block {block}: vertex count");
        assert_eq!(sparse.triangles.len(), dense.triangles.len(), "block {block}: triangle count");
        assert_eq!(sparse.triangles, dense.triangles, "block {block}: triangle indices");
        for (a, b) in sparse.vertices.iter().zip(&dense.vertices) {
            assert_eq!(a, b, "block {block}: vertex position");
        }
        println!(
            "block {block}: kept {}/{} blocks, {} of {} evaluations ({:.1} %)",
            stats.blocks_kept,
            stats.blocks,
            stats.evaluations,
            stats.dense_evaluations,
            100.0 * stats.work_fraction()
        );
    }
}

#[test]
fn culling_work_falls_toward_the_surface_area() {
    // A sphere's boundary is a thin shell of its box, so the saving grows
    // with the resolution: the dense cost is the volume, the sparse cost
    // the area.
    let sphere = Solid::sphere([0.0; 3], 0.5);
    let mut fractions = Vec::new();
    for div in [16usize, 32, 64] {
        let grid = Grid3::covering(sphere.bbox(), 0.5 / div as f64, 2).unwrap();
        let (mesh, stats) = surface_nets_sparse(&sphere, &grid, 2, 8);
        assert!(mesh.report().watertight);
        println!(
            "r/h = {div}: {:.1} % of the dense evaluations, {} of {} blocks kept",
            100.0 * stats.work_fraction(),
            stats.blocks_kept,
            stats.blocks
        );
        fractions.push(stats.work_fraction());
    }
    assert!(
        fractions[2] < fractions[0],
        "the saving should grow with resolution: {fractions:?}"
    );
    assert!(fractions[2] < 0.6, "finest work fraction {:.3}", fractions[2]);

    // And the mesh is still the exact sphere to the usual order.
    let grid = Grid3::covering(sphere.bbox(), 0.5 / 32.0, 2).unwrap();
    let (mesh, _) = surface_nets_sparse(&sphere, &grid, 2, 8);
    let exact = 4.0 / 3.0 * PI * 0.5f64.powi(3);
    assert!((mesh.report().volume - exact).abs() / exact < 5e-4);
}

#[test]
fn anisotropic_cells_resolve_a_thin_feature_in_a_wide_model() {
    // A plate 1000 times wider than it is thick. An isotropic grid fine
    // enough for the thickness needs 1e9 cells; 2 um in z and 200 um in x
    // and y needs 1e5 and resolves it exactly as well.
    let t = 20e-6;
    let plate = Solid::box_between([-0.01, -0.01, 0.0], [0.01, 0.01, t]);
    let grid = Grid3::covering_with(plate.bbox(), [200e-6, 200e-6, 4e-6], 2).unwrap();
    println!(
        "anisotropic grid: {} x {} x {} = {} nodes",
        grid.dims[0], grid.dims[1], grid.dims[2], grid.nodes()
    );
    assert!(grid.nodes() < 200_000, "the grid should be small: {}", grid.nodes());
    let mesh = surface_nets(&plate, &grid, 2);
    let report = mesh.report();
    assert!(report.watertight, "{report:?}");
    let exact = 0.02 * 0.02 * t;
    let rel = (report.volume - exact).abs() / exact;
    println!("plate volume {:e} vs exact {exact:e} ({rel:.3e})", report.volume);
    assert!(rel < 0.02, "volume off by {rel:e}");
    // The thin direction really is resolved: the mesh spans the thickness.
    let (lo, hi) = report.bbox;
    assert!(((hi[2] - lo[2]) - t).abs() / t < 0.02, "z extent {:e} vs {t:e}", hi[2] - lo[2]);

    // The same cell count isotropically cannot resolve it at all: at
    // 200 um the 20 um plate is thinner than one cell and the surface is
    // lost.
    let coarse = mesh_solid(&plate, 200e-6).unwrap();
    assert!(
        coarse.triangles.is_empty() || coarse.report().volume < 0.5 * exact,
        "an isotropic cell coarser than the feature must not find it"
    );

    assert!(mesh_solid_anisotropic(&plate, [200e-6, 200e-6, 4e-6]).is_ok());
    assert!(mesh_solid_anisotropic(&plate, [200e-6, 200e-6, 0.0]).is_err());
}

#[test]
fn a_whole_pcb_meshes_with_the_two_together() {
    // The case `analysis/cad.md` listed as out of reach: 35 um copper on a
    // 1.5 mm board, 40 mm square. Isotropic cells fine enough for the
    // copper would need about 1e10 nodes. Anisotropic cells bring that to
    // 1e7, and block culling evaluates a fraction of those.
    let half = 0.02;
    let board = Board::new(
        Region2::rect([0.0, 0.0], [half, half]),
        vec![
            PcbLayer::Copper { ounces: 1.0 },
            PcbLayer::Dielectric {
                material: "FR-4".into(),
                thickness: 1.5e-3,
                eps_r: 4.4,
                loss_tangent: 0.02,
                kappa: 0.29,
            },
            PcbLayer::Copper { ounces: 1.0 },
        ],
    )
    .with(Feature::Trace {
        layer: 1,
        points: vec![[-0.015, 0.0], [0.0, 0.0], [0.0, 0.012], [0.015, 0.012]],
        width: 0.6e-3,
    })
    .with(Feature::Pour { layer: 0, region: Region2::rect([0.0, 0.0], [half, half]) });

    let solid = board.solid();
    let t = board.thickness();
    // 12 um through the thickness resolves the 34.8 um copper; 250 um in
    // the plane resolves a 600 um trace.
    let spacing = [250e-6, 250e-6, 12e-6];
    let grid = Grid3::covering_with(solid.bbox(), spacing, 2).unwrap();
    let isotropic_nodes = (0.04 / 12e-6_f64).powi(2) * (t / 12e-6);
    println!(
        "board grid {} x {} x {} = {} nodes; isotropic at 12 um would need {:.1e}",
        grid.dims[0], grid.dims[1], grid.dims[2], grid.nodes(), isotropic_nodes
    );
    assert!(grid.nodes() < 4_000_000, "{} nodes", grid.nodes());
    assert!(isotropic_nodes / grid.nodes() as f64 > 100.0, "the saving should be large");

    let (mesh, stats) = surface_nets_sparse(&solid, &grid, 2, 8);
    println!(
        "  meshed: {} triangles, {} of {} blocks kept, {:.1} % of the evaluations",
        mesh.triangles.len(),
        stats.blocks_kept,
        stats.blocks,
        100.0 * stats.work_fraction()
    );
    // Culling buys nothing here, and that is the expected answer: a board
    // is a slab, so its surfaces pass through every block in the plane.
    // The two techniques address different shapes — anisotropy a thin
    // feature, culling a compact one — and this case is won entirely by
    // the anisotropy.
    assert!(
        stats.work_fraction() > 0.9,
        "a slab leaves nothing to cull: {:.3}",
        stats.work_fraction()
    );
    let report = mesh.report();
    assert_eq!(report.boundary_edges, 0, "the board surface must be closed: {report:?}");
    assert!(report.volume > 0.0);
    // Dielectric plus two copper layers, where the copper exists.
    let dielectric = 0.04 * 0.04 * 1.5e-3;
    assert!(
        report.volume > dielectric && report.volume < dielectric * 1.1,
        "volume {:e} against the dielectric's {dielectric:e}",
        report.volume
    );
    // The copper layers are genuinely in the mesh: the surface spans the
    // whole stack, not just the dielectric.
    let (lo, hi) = report.bbox;
    assert!(
        ((hi[2] - lo[2]) - t).abs() < 2.0 * 12e-6,
        "z extent {:e} vs the stack's {t:e}",
        hi[2] - lo[2]
    );
    assert!(OZ_COPPER > 12e-6, "the copper must be thicker than one cell to be resolved");
}

#[test]
fn the_grid_reports_what_it_can_resolve() {
    let g = Grid3::covering_with(([0.0; 3], [1.0, 2.0, 3.0]), [0.1, 0.2, 0.05], 0).unwrap();
    assert!((g.finest() - 0.05).abs() < 1e-18);
    let expect = 0.5 * (0.1f64.powi(2) + 0.2f64.powi(2) + 0.05f64.powi(2)).sqrt();
    assert!((g.half_diagonal() - expect).abs() < 1e-18);
    assert_eq!(g.cells(), [g.dims[0] - 1, g.dims[1] - 1, g.dims[2] - 1]);
    let (lo, hi) = g.bounds();
    assert!(lo[0].abs() < 1e-18 && hi[0] >= 1.0 - 1e-12);
    assert!(Grid3::covering_with(([0.0; 3], [1.0; 3]), [0.1, -0.1, 0.1], 0).is_err());
}
