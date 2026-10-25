//! The vapour-compression cycle against the same cycle computed by
//! CoolProp (5 °C evaporating, 45 °C condensing, 5 K superheat, 3 K
//! subcooling, isentropic efficiency 0.7).

use hvac::refrigeration::{cycle, refrigerant};

// (fluid, p_evap, p_cond, h1, h2, h3, T2) from CoolProp 8.0.0 PropsSI.
const COOLPROP: [(&str, f64, f64, f64, f64, f64, f64); 5] = [
    ("R134a", 349658.6078613138, 1159924.238342344, 406070.70396711427, 442671.877591148, 259386.06281552024, 336.90315078452875),
    ("R32", 951448.0196910894, 2794781.0255771917, 522425.54869951465, 586249.0957687902, 279608.18765672995, 368.16287466063307),
    ("R290", 551116.6057802361, 1534314.068491049, 589259.9610116436, 659525.5569962905, 312893.0425635099, 336.32225142139123),
    ("R1234yf", 373006.0804758775, 1153889.3384393516, 371447.3795297649, 400992.0710431075, 257961.77826324612, 327.53279534467083),
    ("R410A", 933176.1684178826, 2733757.479553493, 428538.48883734277, 471651.1131727984, 269798.76566395117, 350.1484799152381),
];

#[test]
fn cycles_match_coolprop() {
    let (te, tc) = (278.15, 318.15);
    for (name, pe, pc, h1, h2, h3, t2) in COOLPROP {
        let eos = refrigerant(name).unwrap();
        let c = cycle(&eos, te, tc, 5.0, 3.0, 0.7).unwrap();
        let cop_ref = (h1 - h3) / (h2 - h1);
        println!(
            "{name}: COP {:.5} (CoolProp {cop_ref:.5}), p_evap {:.1} kPa ({:.1}), p_cond {:.1} ({:.1}), discharge {:.2} °C ({:.2}), inlet quality {:.3}",
            c.cop_cooling, c.p_evap / 1e3, pe / 1e3, c.p_cond / 1e3, pc / 1e3, c.t_discharge - 273.15, t2 - 273.15, c.quality_in
        );
        // Carnot bounds the cooling COP from above.
        assert!(c.cop_cooling < te / (tc - te));
        assert!((c.cop_heating - c.cop_cooling - 1.0).abs() < 1e-12);
        if name == "R410A" {
            // A blend: CoolProp's dew (evaporator) and bubble (condenser)
            // pressures come from ancillary fits; the pseudo-pure equation
            // has one saturation pressure between them, so the cycle
            // agrees to its glide, not to round-off.
            assert!((c.cop_cooling / cop_ref - 1.0).abs() < 5e-3);
        } else {
            assert!((c.p_evap / pe - 1.0).abs() < 1e-9 && (c.p_cond / pc - 1.0).abs() < 1e-9);
            assert!((c.h[0] - h1).abs() < 1e-3 && (c.h[1] - h2).abs() < 5e-3 && (c.h[2] - h3).abs() < 1e-3);
            assert!((c.t_discharge - t2).abs() < 1e-6);
            println!("   dh1 {:.2e} dh2 {:.2e} dh3 {:.2e} J/kg, dT2 {:.2e} K, COP rel {:.2e}", c.h[0] - h1, c.h[1] - h2, c.h[2] - h3, c.t_discharge - t2, c.cop_cooling / cop_ref - 1.0);
            // COP divides by the compressor work h₂ − h₁, a small difference
            // of enthalpies that agree to 1e-9 (the property test's bound):
            // the COP bound is that, amplified by h₂/(h₂ − h₁).
            let amplification = h2 / (h2 - h1);
            assert!((c.cop_cooling / cop_ref - 1.0).abs() < 1e-9 * amplification, "COP");
        }
    }
}

#[test]
fn refusals() {
    let co2 = refrigerant("R744").unwrap();
    let e = cycle(&co2, 273.15, 313.15, 5.0, 0.0, 0.7).unwrap_err();
    assert!(e.contains("transcritical"), "{e}");
    assert!(refrigerant("R22").is_err());
}
