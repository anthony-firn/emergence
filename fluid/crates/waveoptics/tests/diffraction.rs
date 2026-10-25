//! The diffraction ladder against the exact Rayleigh–Sommerfeld integral.

use num_complex::Complex64 as Cx;
use std::f64::consts::PI;
use waveoptics::diffraction::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    }
}

const LAMBDA: f64 = 632.8e-9; // He–Ne

fn random_query(rng: &mut Rng) -> Query {
    let ax = LAMBDA * (12.0 + 30.0 * rng.next());
    let ay = LAMBDA * (12.0 + 30.0 * rng.next());
    let l = 10f64.powf(-4.0 + 2.5 * rng.next()); // 0.1 mm … 32 mm
    let incident = if rng.next() < 0.5 {
        Incident::Plane {
            amplitude: Cx::new(1.0, 0.0),
            sines: [0.05 * (rng.next() - 0.5), 0.05 * (rng.next() - 0.5)],
        }
    } else {
        Incident::Point {
            amplitude: Cx::new(1e-3, 0.0),
            position: [1e-4 * (rng.next() - 0.5), 1e-4 * (rng.next() - 0.5)],
            distance: 10f64.powf(-3.0 + 2.0 * rng.next()),
        }
    };
    // Observation points up to ~0.4 rad off axis.
    let off = 0.4 * l;
    Query {
        aperture: Rect {
            center: [2e-5 * (rng.next() - 0.5), 2e-5 * (rng.next() - 0.5)],
            size: [ax, ay],
        },
        incident,
        point: [off * (2.0 * rng.next() - 1.0), off * (2.0 * rng.next() - 1.0), l],
        wavelength: LAMBDA,
    }
}

#[test]
fn rigorous_claims_hold_against_rayleigh_sommerfeld() {
    use rayon::prelude::*;
    let mut rng = Rng(7);
    let qs: Vec<Query> = (0..300).map(|_| random_query(&mut rng)).collect();
    let rows: Vec<(f64, f64, f64, f64, f64)> = qs
        .par_iter()
        .map(|q| {
            let r3 = rayleigh_sommerfeld(q, 3);
            let r6 = rayleigh_sommerfeld(q, 6);
            let ref_err = (r3 - r6).norm();
            let (fq, fl) = (local_expansion(q, Order::Quadratic), local_expansion(q, Order::Linear));
            // the classic centre-corrected Fresnel and its bound, too
            let (bc, _) = bounds(q);
            assert!((fresnel(q) - r6).norm() <= bc + ref_err, "classic Fresnel bound violated");
            let li = line_integral(q);
            assert!((li.u - r6).norm() <= li.bound + ref_err + 1e-13 * r6.norm(), "line-integral bound violated: {} > {}", (li.u - r6).norm(), li.bound);
            ((fq.u - r6).norm(), fq.bound, (fl.u - r6).norm(), fl.bound, ref_err)
        })
        .collect();
    let (mut wf, mut wfr, mut wref) = (0.0f64, 0.0f64, 0.0f64);
    for (i, &(ef, bf, efr, bfr, re)) in rows.iter().enumerate() {
        assert!(ef <= bf + re, "query {i}: Fresnel error {ef:e} > bound {bf:e} (+ref {re:e}) {:?}", qs[i]);
        assert!(efr <= bfr + re, "query {i}: Fraunhofer error {efr:e} > bound {bfr:e} (+ref {re:e}) {:?}", qs[i]);
        wf = wf.max(ef / bf);
        wfr = wfr.max(efr / bfr);
        wref = wref.max(re);
    }
    println!("300 queries: worst error/bound tilted Fresnel {wf:.3}, tilted Fraunhofer {wfr:.3}; worst reference uncertainty {wref:.1e}");
}

#[test]
fn far_field_matches_the_analytic_sinc_and_conserves_power() {
    // Normal plane wave on a 50 µm × 40 µm opening, screen at 0.5 m.
    let ap = Rect { center: [0.0; 2], size: [50e-6, 40e-6] };
    let l = 0.5;
    let q = |x: f64, y: f64| Query {
        aperture: ap,
        incident: Incident::Plane { amplitude: Cx::new(1.0, 0.0), sines: [0.0; 2] },
        point: [x, y, l],
        wavelength: LAMBDA,
    };
    // Analytic Fraunhofer: |U| = (A a b/(λ r)) (L/r) |sinc(k a x/(2r)) sinc(k b y/(2r))|.
    let k = 2.0 * PI / LAMBDA;
    for &(x, y) in &[(0.0, 0.0), (3e-3, 0.0), (6.3e-3, 2e-3), (-1e-2, 5e-3)] {
        let r = (x * x + y * y + l * l).sqrt();
        let s = |v: f64| if v == 0.0 { 1.0 } else { v.sin() / v };
        let analytic = ap.area() / (LAMBDA * r) * (l / r) * (s(k * 50e-6 * x / (2.0 * r)) * s(k * 40e-6 * y / (2.0 * r))).abs();
        let u = fraunhofer(&q(x, y)).norm();
        assert!((u - analytic).abs() <= 1e-6 * analytic.max(1e-3), "({x}, {y}): {u} vs {analytic}");
    }
    // Power through the opening (|A|² × area) equals the far-field power
    // (integrated over a screen out to the 40th zero in x and y).
    let (half_x, half_y) = (40.0 * LAMBDA * l / 50e-6, 40.0 * LAMBDA * l / 40e-6);
    let n = 800;
    let mut p = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x = -half_x + 2.0 * half_x * (i as f64 + 0.5) / n as f64;
            let y = -half_y + 2.0 * half_y * (j as f64 + 0.5) / n as f64;
            let qq = q(x, y);
            p += fresnel(&qq).norm_sqr() * (2.0 * half_x / n as f64) * (2.0 * half_y / n as f64);
        }
    }
    let rel = (p - ap.area()) / ap.area();
    println!("far-field power / aperture power − 1 = {rel:.2e}");
    // The sinc² tails beyond the 40th zero carry ~1/(π² 40)·2 per axis ≈ 1 %.
    assert!(rel < 0.0 && rel > -0.015, "{rel}");
}

#[test]
fn ladder_uses_cheap_rungs_where_certified() {
    let ladder = ladder(0);
    // Lab regime: 50 µm × 3 mm slits, screens 5 cm … 2 m away, points over
    // a 4 cm field. Near and far, the ladder must certify cheap rungs.
    let mut rng = Rng(3);
    let qs: Vec<Query> = (0..200)
        .map(|_| {
            let l = 0.05 * 40f64.powf(rng.next());
            Query {
                aperture: Rect { center: [0.0, 0.0], size: [50e-6, 3e-3] },
                incident: Incident::Plane { amplitude: Cx::new(1.0, 0.0), sines: [0.0; 2] },
                point: [0.02 * (2.0 * rng.next() - 1.0), 0.004 * (2.0 * rng.next() - 1.0), l],
                wavelength: LAMBDA,
            }
        })
        .collect();
    let scale = |q: &Query| {
        let a = match q.incident {
            Incident::Plane { amplitude, .. } => amplitude.norm(),
            Incident::Point { amplitude, distance, .. } => amplitude.norm() / distance,
        };
        a * q.aperture.area() / (q.wavelength * q.point[2])
    };
    let mut counts = [0usize; 4];
    for q in &qs {
        let s = ladder.solve(q, emergence::Request::new(1e-3 * scale(q)));
        assert!(s.met);
        counts[s.rung] += 1;
    }
    println!("3 mm slit, 5 cm … 2 m: rungs (fraunhofer, fresnel, line, rs) {counts:?}");
    // The 3 mm slit is ~3.5 Fresnel zones tall at 1 m: Fraunhofer cannot be
    // certified; Fresnel or the line integral must cover every point.
    assert_eq!(counts[0], 0);
    assert_eq!(counts[3], 0, "no point should need the full 2D integral");
    // A 20 µm square 2 m away is in the far field.
    let far = Query {
        aperture: Rect { center: [0.0; 2], size: [20e-6, 20e-6] },
        incident: Incident::Plane { amplitude: Cx::new(1.0, 0.0), sines: [0.0; 2] },
        point: [0.01, 0.005, 2.0],
        wavelength: LAMBDA,
    };
    let s = ladder.solve(&far, emergence::Request::new(1e-3 * scale(&far)));
    assert_eq!(ladder.rung_name(s.rung), "tilted fraunhofer");
}
