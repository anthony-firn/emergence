//! Geometry checks against analytic shapes, brute force, and the reference
//! VTK library.

use std::f64::consts::PI;

use geometry::bvh::{closest_on_triangle, Bvh};
use geometry::primitives::*;
use geometry::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    }
    fn point(&mut self, lo: f64, hi: f64) -> P3 {
        [lo + (hi - lo) * self.next(), lo + (hi - lo) * self.next(), lo + (hi - lo) * self.next()]
    }
}

fn tmp(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("geometry-test-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d.join(name)
}

#[test]
fn primitives_are_watertight_with_converging_mass_properties() {
    let cube = cuboid([0.0; 3], [1.0, 2.0, 3.0]);
    let r = cube.report();
    assert!(r.watertight);
    assert!((r.volume - 6.0).abs() < 1e-14 && (r.area - 22.0).abs() < 1e-13);
    let c = cube.centroid();
    assert!((c[0] - 0.5).abs() < 1e-14 && (c[1] - 1.0).abs() < 1e-14 && (c[2] - 1.5).abs() < 1e-14);

    // Icosphere: polyhedral volume converges to 4π/3 at O(h²) (×4 per level).
    let mut prev = f64::INFINITY;
    for level in 2..6 {
        let s = icosphere([0.0; 3], 1.0, level);
        let r = s.report();
        assert!(r.watertight, "level {level}");
        let err = (r.volume - 4.0 * PI / 3.0).abs() / (4.0 * PI / 3.0);
        assert!(err < prev / 3.5, "level {level}: {err:e} vs {prev:e}");
        prev = err;
    }
    assert!(prev < 1e-3);

    // Torus: V = 2π² R r², A = 4π² R r (limits).
    let t = torus([0.0; 3], 2.0, 0.5, 256, 128);
    let r = t.report();
    assert!(r.watertight);
    assert!((r.volume - 2.0 * PI * PI * 2.0 * 0.25).abs() / (PI * PI) < 1e-3);
    assert!((r.area - 4.0 * PI * PI).abs() / (4.0 * PI * PI) < 1e-3);
    // Cylinder.
    let c = cylinder_z([0.0; 3], 0.5, 2.0, 512);
    assert!(c.report().watertight);
    assert!((c.volume() - PI * 0.25 * 2.0).abs() / (PI * 0.5) < 1e-4);
}

#[test]
fn stl_and_obj_round_trip() {
    let s = icosphere([0.1, -0.2, 0.3], 0.7, 3);
    let path = tmp("sphere.stl");
    s.write_stl(&path).unwrap();
    let back = TriMesh::read(&path).unwrap();
    let (a, b) = (s.report(), back.report());
    assert!(b.watertight);
    assert_eq!((a.vertices, a.triangles), (b.vertices, b.triangles));
    // f32 storage: relative volume change ~1e-7.
    assert!((a.volume - b.volume).abs() / a.volume < 1e-6);

    let ascii = "solid t\nfacet normal 0 0 1\n outer loop\n  vertex 0 0 0\n  vertex 1 0 0\n  vertex 0 1 0\n endloop\nendfacet\nendsolid t\n";
    let m = TriMesh::from_stl_bytes(ascii.as_bytes()).unwrap();
    assert_eq!(m.triangles.len(), 1);
    assert!((m.area() - 0.5).abs() < 1e-15);

    // OBJ cube with quads and v/vt/vn indices.
    let obj = "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 0 0 1\nv 1 0 1\nv 1 1 1\nv 0 1 1\n\
               f 1/1/1 4 3 2\nf 5 6 7 8\nf 1 2 6 5\nf 2 3 7 6\nf 3 4 8 7\nf 4 1 5 8\n";
    let mut m = TriMesh::from_obj(obj).unwrap();
    let r = m.report();
    assert!(r.watertight && r.triangles == 12);
    m.orient_outward();
    assert!((m.volume() - 1.0).abs() < 1e-15);
}

/// Exact signed distance to the axis-aligned box `[lo, hi]`.
fn box_sdf(p: P3, lo: P3, hi: P3) -> f64 {
    let c = [(lo[0] + hi[0]) / 2.0, (lo[1] + hi[1]) / 2.0, (lo[2] + hi[2]) / 2.0];
    let e = [(hi[0] - lo[0]) / 2.0, (hi[1] - lo[1]) / 2.0, (hi[2] - lo[2]) / 2.0];
    let q = [(p[0] - c[0]).abs() - e[0], (p[1] - c[1]).abs() - e[1], (p[2] - c[2]).abs() - e[2]];
    let outside = norm([q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)]);
    outside + q[0].max(q[1]).max(q[2]).min(0.0)
}

#[test]
fn box_signed_distance_is_exact_including_edges_and_corners() {
    let (lo, hi) = ([-0.5, -1.0, -0.25], [0.5, 1.0, 0.25]);
    let sdf = Sdf::new(cuboid(lo, hi));
    assert_eq!(sdf.sign, SignMode::PseudoNormal);
    let mut rng = Rng(1);
    for _ in 0..20000 {
        let p = rng.point(-2.0, 2.0);
        let d = sdf.distance(p);
        assert!((d - box_sdf(p, lo, hi)).abs() < 1e-12, "p {p:?}: {d} vs {}", box_sdf(p, lo, hi));
    }
    // Points whose closest feature is an edge or a corner (the hard cases
    // for the sign): just outside, diagonal from the corners.
    for s in [[1.0, 1.0, 1.0], [-1.0, 1.0, -1.0], [1.0, -1.0, 0.0], [0.0, 1.0, 1.0]] {
        let p = [s[0] * 0.5 + 1e-6 * s[0], s[1] * 1.0 + 1e-6 * s[1], s[2] * 0.25 + 1e-6 * s[2]];
        assert!(sdf.distance(p) > 0.0, "{p:?}");
    }
}

#[test]
fn bvh_closest_point_matches_brute_force() {
    let mesh = torus([0.3, 0.0, -0.1], 1.0, 0.3, 48, 24);
    let bvh = Bvh::build(&mesh);
    let mut rng = Rng(5);
    for _ in 0..3000 {
        let p = rng.point(-2.0, 2.0);
        let hit = bvh.closest(&mesh, p).unwrap();
        let brute = (0..mesh.triangles.len())
            .map(|t| {
                let [a, b, c] = mesh.corners(t);
                let q = closest_on_triangle(p, a, b, c).0;
                dot(sub(p, q), sub(p, q))
            })
            .fold(f64::INFINITY, f64::min);
        assert!((hit.dist2 - brute).abs() <= 1e-14 * (1.0 + brute));
    }
}

#[test]
fn curved_surfaces_within_their_tessellation_error() {
    // Sphere: the polyhedron lies between r cos(θ_max) and r.
    let r = 0.8;
    let sphere = icosphere([0.0; 3], r, 5);
    // Hausdorff distance between the inscribed polyhedron and the sphere:
    // the largest gap r − dist(centre, face) over all faces.
    let sagitta = (0..sphere.triangles.len())
        .map(|t| {
            let [a, b, c] = sphere.corners(t);
            let q = closest_on_triangle([0.0; 3], a, b, c).0;
            r - norm(q)
        })
        .fold(0.0, f64::max);
    let sdf = Sdf::new(sphere);
    let mut rng = Rng(9);
    let mut worst: f64 = 0.0;
    for _ in 0..20000 {
        let p = rng.point(-1.5, 1.5);
        let exact = norm(p) - r;
        let d = sdf.distance(p);
        worst = worst.max((d - exact).abs());
        if exact.abs() > sagitta {
            assert_eq!(d < 0.0, exact < 0.0, "sign at {p:?}");
        }
    }
    assert!(worst <= sagitta * (1.0 + 1e-9), "sphere worst {worst:e} vs Hausdorff {sagitta:e}");
    // Torus: analytic SDF √((ρ − R)² + z²) − r.
    let (big, small) = (1.0, 0.35);
    let sdf = Sdf::new(torus([0.0; 3], big, small, 400, 160));
    let mut worst: f64 = 0.0;
    for _ in 0..20000 {
        let p = rng.point(-1.6, 1.6);
        let rho = (p[0] * p[0] + p[1] * p[1]).sqrt();
        let exact = ((rho - big).powi(2) + p[2] * p[2]).sqrt() - small;
        worst = worst.max((sdf.distance(p) - exact).abs());
    }
    // Chord sagitta of 400 × 160 segments: (R + r)(1 − cos(π/400)) + r(1 − cos(π/160)).
    let bound = (big + small) * (1.0 - (PI / 400.0).cos()) + small * (1.0 - (PI / 160.0).cos());
    assert!(worst <= bound * 1.01, "torus worst {worst:e} vs {bound:e}");
}

#[test]
fn open_surface_falls_back_to_winding_numbers() {
    let mut s = icosphere([0.0; 3], 1.0, 3);
    s.triangles.truncate(s.triangles.len() - 6); // punch a hole
    let sdf = Sdf::new(s);
    assert_eq!(sdf.sign, SignMode::WindingNumber);
    let mut rng = Rng(3);
    for _ in 0..2000 {
        let p = rng.point(-1.4, 1.4);
        let exact = norm(p) - 1.0;
        if exact.abs() > 0.05 {
            assert_eq!(sdf.distance(p) < 0.0, exact < 0.0, "{p:?}");
        }
    }
}

/// Python with the `vtk` module, if available.
fn python_vtk() -> bool {
    std::process::Command::new("python3")
        .args(["-c", "import vtk"])
        .output()
        .is_ok_and(|o| o.status.success())
}

#[test]
fn vtk_files_read_back_with_the_reference_library() {
    if !python_vtk() {
        eprintln!("python3 with vtk not available: skipping");
        return;
    }
    let cells = [4, 3, 2];
    let n: usize = cells.iter().product();
    let pressure: Vec<f64> = (0..n).map(|i| i as f64 * 0.5 - 1.0).collect();
    let vel: Vec<P3> = (0..n).map(|i| [i as f64, -(i as f64), 0.25]).collect();
    let np: usize = cells.iter().map(|c| c + 1).product();
    let sdfv: Vec<f64> = (0..np).map(|i| (i as f64).sqrt()).collect();
    let grid = tmp("grid.vtk");
    vtk::write_grid(&grid, [0.1, 0.2, 0.3], [0.5, 0.25, 1.0], cells, &[("p", vtk::Data::Scalars(&pressure)), ("u", vtk::Data::Vectors(&vel))], &[("phi", vtk::Data::Scalars(&sdfv))]).unwrap();
    let mesh = icosphere([0.0; 3], 1.0, 1);
    let cp: Vec<f64> = (0..mesh.triangles.len()).map(|t| t as f64).collect();
    let pn: Vec<P3> = mesh.vertices.clone();
    let surf = tmp("surf.vtk");
    vtk::write_surface(&surf, &mesh, &[("n", vtk::Data::Vectors(&pn))], &[("cp", vtk::Data::Scalars(&cp))]).unwrap();
    let script = format!(
        r#"
import vtk
from vtk.util.numpy_support import vtk_to_numpy as v
r = vtk.vtkStructuredPointsReader(); r.SetFileName("{g}"); r.Update(); o = r.GetOutput()
assert o.GetDimensions() == (5, 4, 3), o.GetDimensions()
assert abs(o.GetSpacing()[1] - 0.25) < 1e-15 and abs(o.GetOrigin()[2] - 0.3) < 1e-15
p = v(o.GetCellData().GetArray("p")); u = v(o.GetCellData().GetArray("u")); f = v(o.GetPointData().GetArray("phi"))
assert p.shape == (24,) and abs(p[5] - 1.5) < 1e-15, p
assert u.shape == (24, 3) and u[7][1] == -7.0 and u[3][2] == 0.25
assert abs(f[59] - 59 ** 0.5) < 1e-15
s = vtk.vtkPolyDataReader(); s.SetFileName("{s}"); s.Update(); m = s.GetOutput()
assert m.GetNumberOfPoints() == {nv} and m.GetNumberOfCells() == {nt}
assert v(m.GetCellData().GetArray("cp"))[13] == 13.0
assert abs(v(m.GetPointData().GetArray("n"))[2][0] - ({x})) < 1e-15
print("ok")
"#,
        g = grid.display(),
        s = surf.display(),
        nv = mesh.vertices.len(),
        nt = mesh.triangles.len(),
        x = mesh.vertices[2][0]
    );
    let out = std::process::Command::new("python3").args(["-c", &script]).output().unwrap();
    assert!(out.status.success(), "vtk read-back failed: {}", String::from_utf8_lossy(&out.stderr));
}

#[test]
fn step_import_through_gmsh_matches_the_analytic_solid() {
    let Some(python) = geometry::cad::gmsh_python() else {
        eprintln!("gmsh not set up (fluid/scripts/setup-cad.sh): skipping");
        return;
    };
    // Make a STEP cylinder r = 0.05, L = 0.41 with OpenCASCADE, then import.
    let step = tmp("cylinder.step");
    let make = format!(
        "import gmsh\ngmsh.initialize()\ngmsh.option.setNumber('General.Terminal', 0)\n\
         gmsh.model.occ.addCylinder(0.2, 0.2, 0.0, 0, 0, 0.41, 0.05)\ngmsh.model.occ.synchronize()\n\
         gmsh.write('{}')\ngmsh.finalize()\n",
        step.display()
    );
    let out = std::process::Command::new(&python).args(["-c", &make]).output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let mesh = TriMesh::read(&step).unwrap();
    let r = mesh.report();
    assert!(r.watertight, "{r:?}");
    let exact = PI * 0.05 * 0.05 * 0.41;
    let rel = (r.volume.abs() - exact).abs() / exact;
    // Inscribed polygon with n = 64 per 2π: area deficit 1 − sin(2π/n) n/(2π) ≈ 1.6e-3.
    let deficit = 1.0 - (2.0 * PI / 64.0).sin() * 64.0 / (2.0 * PI);
    println!("STEP cylinder: {} triangles, volume error {rel:.2e} (inscribed-polygon deficit {deficit:.2e})", r.triangles);
    assert!(rel <= 1.2 * deficit, "volume error {rel:e}");
    let (lo, hi) = r.bbox;
    assert!((lo[2]).abs() < 1e-6 && (hi[2] - 0.41).abs() < 1e-6 && (hi[0] - 0.25).abs() < 1e-6);
}
