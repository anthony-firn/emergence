//! The JSON case runner end to end: geometry file in, forces/summary/VTK out.

use geometry::primitives::cylinder_z;
use immersed::case::Case;

#[test]
fn case_file_runs_end_to_end_and_writes_readable_outputs() {
    let dir = std::env::temp_dir().join(format!("immersed-case-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    cylinder_z([0.2, 0.2, -0.1], 0.05, 0.2, 256).write_stl(dir.join("cyl.stl")).unwrap();
    let json = r#"{
      "name": "coarse channel",
      "fluid": { "density": 1.2, "kinematic_viscosity": 0.001 },
      "domain": { "min": [0, 0, 0], "max": [1.0, 0.41, 0.01], "cell": 0.01 },
      "boundaries": {
        "x_min": { "type": "inlet", "velocity": [0.3, 0, 0], "profile": "parabolic" },
        "x_max": { "type": "outlet" },
        "y_min": { "type": "wall" }, "y_max": { "type": "wall" },
        "z_min": { "type": "periodic" }, "z_max": { "type": "periodic" }
      },
      "bodies": [ { "file": "cyl.stl" } ],
      "run": { "t_end": 1.0, "tolerance": 1e-9 },
      "output": { "directory": "out", "reference": { "velocity": 0.2, "area": 0.001 },
                  "probes": [[0.15, 0.2, 0.005]], "every": 0.5 }
    }"#;
    let case = Case::from_json(json).unwrap();
    let s = case.run(&dir, |_| {}).unwrap();
    assert_eq!(s.cells, [100, 41, 1]);
    assert!(!s.converged && (s.time - 1.0).abs() < 1e-9);
    let b = &s.bodies[0];
    assert!(b.watertight && b.force[0] > 0.0, "drag must push downstream: {:?}", b.force);
    // ρ = 1.2 enters the force; the coefficient divides it back out.
    let c = b.coefficients.unwrap();
    assert!((c[0] - b.force[0] / (0.5 * 1.2 * 0.04 * 0.001)).abs() < 1e-12);
    assert!(s.probe_pressure[0] > 0.0, "stagnation pressure");
    let out = dir.join("out");
    let csv = std::fs::read_to_string(out.join("forces.csv")).unwrap();
    assert_eq!(csv.lines().count(), 3, "{csv}");
    let summary: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(out.join("summary.json")).unwrap()).unwrap();
    assert_eq!(summary["cells"][0], 100);
    for f in ["fields.vtk", "body0.vtk"] {
        assert!(out.join(f).metadata().unwrap().len() > 1000, "{f}");
    }
    // Errors are reported, not panics.
    let bad = json.replace("\"max\": [1.0, 0.41, 0.01]", "\"max\": [1.0, 0.413, 0.01]");
    let err = Case::from_json(&bad).unwrap().run(&dir, |_| {}).unwrap_err();
    assert!(err.to_string().contains("whole number of cells"), "{err}");
    let unpaired = json.replace("\"z_max\": { \"type\": \"periodic\" }", "\"z_max\": { \"type\": \"wall\" }");
    assert!(Case::from_json(&unpaired).unwrap().run(&dir, |_| {}).is_err());
}
