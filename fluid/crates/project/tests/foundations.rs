//! Units, interval expressions and the materials library.

use emergence::Rigor;
use project::expr::Scope;
use project::units::{format_si, parse_quantity, parse_unit, Dim};
use project::{Expr, Library, Val};

#[test]
fn units_convert_exactly_and_refuse_mismatches() {
    let (q, _) = parse_quantity("0.5 mm").unwrap();
    assert_eq!(q.si, 0.5e-3);
    assert_eq!(q.dim, Dim::LENGTH);
    // The international inch is exact.
    assert_eq!(parse_quantity("1 in").unwrap().0.si, 0.0254);
    assert!((parse_quantity("4 mil").unwrap().0.si - 101.6e-6).abs() < 1e-18);
    // Compound units, both spellings.
    let k1 = parse_unit("W/(m*K)").unwrap();
    let k2 = parse_unit("W/m/K").unwrap();
    assert_eq!(k1.dim, k2.dim);
    assert_eq!(k1.dim.symbol(), "W/(K·m)");
    assert_eq!(parse_unit("K/W").unwrap().dim.symbol(), "K/W");
    // Prefixes, and exact symbols winning over prefixed ones.
    assert_eq!(parse_unit("min").unwrap().scale, 60.0);
    assert_eq!(parse_unit("mOhm").unwrap().scale, 1e-3);
    assert_eq!(parse_unit("kg").unwrap().scale, 1.0);
    // Celsius is affine and only valid alone.
    assert!((parse_quantity("25 degC").unwrap().0.si - 298.15).abs() < 1e-12);
    assert!(parse_unit("W/degC").is_err());
    assert!(parse_unit("furlong").is_err());
    assert_eq!(format_si(298.15, Dim::TEMPERATURE, Some("degC")), "25 degC");
}

struct Vars(Vec<(&'static str, Val)>);
impl Scope for Vars {
    fn var(&self, id: &str) -> Result<Val, String> {
        self.0.iter().find(|v| v.0 == id).map(|v| v.1).ok_or_else(|| format!("no {id}"))
    }
    fn table(&self, id: &str, _: &Val) -> Result<Val, String> {
        Err(format!("no table {id}"))
    }
}

#[test]
fn expressions_check_dimensions_and_bound_the_result() {
    let w = Val::new(3e-3, 2.9e-3, 3.1e-3, Dim::LENGTH, Rigor::Rigorous).unwrap();
    let t = Val::exact(35e-6, Dim::LENGTH);
    let s = Vars(vec![("trace.w", w), ("trace.t", t)]);
    let area = Expr::parse("trace.w * trace.t + 0.1[mm]^2 * 0").unwrap().eval(&s).unwrap();
    assert_eq!(area.dim, Dim::LENGTH.powi(2));
    assert!((area.lo - 2.9e-3 * 35e-6).abs() < 1e-20 && (area.hi - 3.1e-3 * 35e-6).abs() < 1e-20);
    // Adding a length to an area is caught.
    let err = Expr::parse("trace.w + trace.w * trace.t").unwrap().eval(&s).unwrap_err();
    assert!(err.contains("cannot add"), "{err}");
    // Celsius literals are refused inside expressions (offset ambiguity).
    assert!(Expr::parse("20[degC] + trace.w").is_err());
    // An interval through a reciprocal stays an enclosure.
    let inv = Expr::parse("1 / trace.w").unwrap().eval(&s).unwrap();
    assert!((inv.lo - 1.0 / 3.1e-3).abs() < 1e-9 && (inv.hi - 1.0 / 2.9e-3).abs() < 1e-9);
    assert_eq!(Expr::parse("mat.air.k(ambient) * x.y").unwrap().deps(), vec!["mat.air.k", "ambient", "x.y"]);
    // x^2 over an interval spanning zero is [0, max], not negative.
    let z = Val::new(0.0, -1.0, 2.0, Dim::NONE, Rigor::Rigorous).unwrap();
    let sq = z.powi(2).unwrap();
    assert_eq!((sq.lo, sq.hi), (0.0, 4.0));
}

#[test]
fn air_table_is_self_consistent() {
    // Incropera Table A.4 lists ν = μ/ρ and Pr = μ c_p / k independently;
    // a transcription error in any column would break these identities.
    let lib = Library::builtin();
    for t in [250.0, 300.0, 350.0, 400.0] {
        let at = Val::exact(t, Dim::TEMPERATURE);
        let g = |p: &str| lib.get("air", p).unwrap().at(&at).unwrap().nom;
        assert!((g("mu") / g("rho") / g("nu") - 1.0).abs() < 2e-3, "nu at {t}");
        assert!((g("mu") * g("cp") / g("k") / g("pr") - 1.0).abs() < 3e-3, "Pr at {t}");
    }
    // Water at 300 K (Table A.6), same identity.
    let w = |p: &str| lib.get("water", p).unwrap().reference().nom;
    assert!((w("mu") * w("cp") / w("k") / w("pr") - 1.0).abs() < 2e-3);
    // Out of table is an error, never an extrapolation.
    assert!(lib.get("air", "k").unwrap().at(&Val::exact(500.0, Dim::TEMPERATURE)).is_err());
    // The interval of a table over a temperature range covers both ends.
    let range = Val::new(320.0, 300.0, 340.0, Dim::TEMPERATURE, Rigor::Rigorous).unwrap();
    let k = lib.get("air", "k").unwrap().at(&range).unwrap();
    assert!(k.lo <= 26.3e-3 && k.hi >= 29.26e-3);
    assert_eq!(k.rigor, Rigor::Validated);
}
