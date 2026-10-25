#!/usr/bin/env python3
"""Generate crates/compressible/src/refrigerant_data.rs from the fluid files
bundled with CoolProp (MIT licence; read locally through the Python package,
CoolProp.get_fluid_param_string(name, "JSON")), with reference states
computed by CoolProp itself for the Rust tests to compare against.

Usage: python3 build_refrigerants.py OUT.rs   (needs `pip install CoolProp`)
"""
import json, sys
import CoolProp
import CoolProp.CoolProp as CP
def f(x):
    return repr(float(x))

def terms(lst):
    out = []
    for t in lst:
        ty = t["type"]
        if ty == "ResidualHelmholtzPower":
            for n, d, tt, l in zip(t["n"], t["d"], t["t"], t["l"]):
                out.append(f"Term::Power {{ n: {f(n)}, d: {f(d)}, t: {f(tt)}, l: {f(l)} }}")
        elif ty == "ResidualHelmholtzGaussian":
            for n, d, tt, eta, eps, beta, gam in zip(t["n"], t["d"], t["t"], t["eta"], t["epsilon"], t["beta"], t["gamma"]):
                out.append(f"Term::Gauss {{ n: {f(n)}, d: {f(d)}, t: {f(tt)}, eta: {f(eta)}, eps: {f(eps)}, beta: {f(beta)}, gamma: {f(gam)} }}")
        elif ty == "ResidualHelmholtzNonAnalytic":
            for n, a, b, bb, c, dd, aa, beta in zip(t["n"], t["a"], t["b"], t["B"], t["C"], t["D"], t["A"], t["beta"]):
                out.append(f"Term::NonAnalytic {{ n: {f(n)}, a: {f(a)}, b: {f(b)}, big_b: {f(bb)}, c: {f(c)}, big_d: {f(dd)}, big_a: {f(aa)}, beta: {f(beta)} }}")
        elif ty == "IdealGasHelmholtzLead":
            out.append(f"Term::Lead {{ a1: {f(t['a1'])}, a2: {f(t['a2'])} }}")
        elif ty == "IdealGasHelmholtzLogTau":
            out.append(f"Term::LogTau {{ a: {f(t['a'])} }}")
        elif ty == "IdealGasHelmholtzPlanckEinstein":
            for n, tt in zip(t["n"], t["t"]):
                out.append(f"Term::PlanckEinstein {{ n: {f(n)}, t: {f(tt)} }}")
        elif ty == "IdealGasHelmholtzPlanckEinsteinFunctionT":
            for n, v in zip(t["n"], t["v"]):
                out.append(f"Term::PlanckEinstein {{ n: {f(n)}, t: {f(v / t['Tcrit'])} }}")
        elif ty == "IdealGasHelmholtzPower":
            for n, tt in zip(t["n"], t["t"]):
                out.append(f"Term::PowerTau {{ n: {f(n)}, t: {f(tt)} }}")
        else:
            raise SystemExit(f"unsupported term {ty}")
    return out

def anc(a):
    return f"Ancillary {{ kind: {a['type']!r}, t_r: {f(a['T_r'])}, reducing: {f(a['reducing_value'])}, tau_r: {str(a.get('using_tau_r', False)).lower()}, n: &[{', '.join(map(f, a['n']))}], t: &[{', '.join(map(f, a['t']))}] }}".replace("'", '"')


FLUIDS = [
    ("R32", "r32", "Tillner-Roth & Yokozeki, J. Phys. Chem. Ref. Data 26, 1273 (1997)."),
    ("R134a", "r134a", "Tillner-Roth & Baehr, J. Phys. Chem. Ref. Data 23, 657 (1994)."),
    ("R410A", "r410a", "Lemmon, Int. J. Thermophys. 24, 991 (2003) (pseudo-pure)."),
    ("R290", "r290", "Lemmon, McLinden & Wagner, J. Chem. Eng. Data 54, 3141 (2009)."),
    ("R1234yf", "r1234yf", "Lemmon & Akasaka, Int. J. Thermophys. 43, 119 (2022)."),
    ("CarbonDioxide", "r744", "Span & Wagner, J. Phys. Chem. Ref. Data 25, 1509 (1996)."),
]

def ideal_terms(lst):
    # CoolProp's reference-state offset a1 + a2·τ is two PowerTau terms.
    out, rest = [], []
    for t in lst:
        if t["type"] == "IdealGasHelmholtzEnthalpyEntropyOffset":
            out.append(f"Term::PowerTau {{ n: {f(t['a1'])}, t: 0.0 }}")
            out.append(f"Term::PowerTau {{ n: {f(t['a2'])}, t: 1.0 }}")
        else:
            rest.append(t)
    return terms(rest) + out

def fluid(name, rust_name, cite):
    d = json.loads(CP.get_fluid_param_string(name, "JSON"))
    if isinstance(d, list):
        d = d[0]
    e = d["EOS"][0]
    m = e["molar_mass"]
    red = e["STATES"]["reducing"]
    a = d["ANCILLARIES"]
    lines = [f"/// {cite}", f"pub fn {rust_name}() -> Helmholtz {{", "    Helmholtz {",
             f'        name: "{name}",',
             f"        molar_mass: {f(m)},", f"        gas_constant: {f(e['gas_constant'])},",
             f"        t_c: {f(red['T'])},", f"        rho_c: {f(red['rhomolar'] * m)},", f"        p_c: {f(red['p'])},",
             f"        t_triple: {f(e['Ttriple'])},",
             "        residual: vec![", *[f"            {t}," for t in terms(e["alphar"])], "        ],",
             "        ideal: vec![", *[f"            {t}," for t in ideal_terms(e["alpha0"])], "        ],",
             f"        anc_p: {anc(a.get('pS', a.get('pL')))},", f"        anc_rho_l: {anc(a['rhoL'])},", f"        anc_rho_v: {anc(a['rhoV'])},",
             "    }", "}", ""]
    return "\n".join(lines)

def reference(name, rust_name):
    """CoolProp's own values: saturation at three temperatures, and
    single-phase (T, rho) states, for the Rust tests."""
    tc = CP.PropsSI("Tcrit", name)
    rows = []
    for frac in (0.70, 0.85, 0.95):
        t = frac * tc
        p = CP.PropsSI("P", "T", t, "Q", 0, name)
        p_dew = CP.PropsSI("P", "T", t, "Q", 1, name)
        rl = CP.PropsSI("Dmass", "T", t, "Q", 0, name)
        rv = CP.PropsSI("Dmass", "T", t, "Q", 1, name)
        hl = CP.PropsSI("Hmass", "T", t, "Q", 0, name)
        hv = CP.PropsSI("Hmass", "T", t, "Q", 1, name)
        sl = CP.PropsSI("Smass", "T", t, "Q", 0, name)
        sv = CP.PropsSI("Smass", "T", t, "Q", 1, name)
        rows.append(f"    Sat {{ t: {f(t)}, p: {f(p)}, p_dew: {f(p_dew)}, rho_l: {f(rl)}, rho_v: {f(rv)}, h_l: {f(hl)}, h_v: {f(hv)}, s_l: {f(sl)}, s_v: {f(sv)} }},")
    states = []
    for (t, p) in ((0.9 * tc, 0.3 * CP.PropsSI("pcrit", name)), (1.2 * tc, 2.0 * CP.PropsSI("pcrit", name)), (0.8 * tc, 0.05 * CP.PropsSI("pcrit", name))):
        rho = CP.PropsSI("Dmass", "T", t, "P", p, name)
        h = CP.PropsSI("Hmass", "T", t, "P", p, name)
        s = CP.PropsSI("Smass", "T", t, "P", p, name)
        cp = CP.PropsSI("Cpmass", "T", t, "P", p, name)
        w = CP.PropsSI("speed_of_sound", "T", t, "P", p, name)
        states.append(f"    State {{ t: {f(t)}, rho: {f(rho)}, p: {f(p)}, h: {f(h)}, s: {f(s)}, cp: {f(cp)}, w: {f(w)} }},")
    up = rust_name.upper()
    return "\n".join([f"pub const {up}_SAT: &[Sat] = &[", *rows, "];", f"pub const {up}_STATES: &[State] = &[", *states, "];", ""])

out = [f"// Generated by fluid/data/build_refrigerants.py from the fluid files bundled with CoolProp {CoolProp.__version__} (MIT). Do not edit.",
       "use crate::helmholtz::{Ancillary, Helmholtz, Term};", ""]
for name, rust, cite in FLUIDS:
    out.append(fluid(name, rust, cite))
out += ["/// CoolProp's values, for the tests: saturation states.",
        "#[derive(Clone, Copy, Debug)]",
        "pub struct Sat { pub t: f64, pub p: f64, pub p_dew: f64, pub rho_l: f64, pub rho_v: f64, pub h_l: f64, pub h_v: f64, pub s_l: f64, pub s_v: f64 }",
        "/// CoolProp's values, for the tests: single-phase states.",
        "#[derive(Clone, Copy, Debug)]",
        "pub struct State { pub t: f64, pub rho: f64, pub p: f64, pub h: f64, pub s: f64, pub cp: f64, pub w: f64 }", ""]
for name, rust, _ in FLUIDS:
    out.append(reference(name, rust))
open(sys.argv[1], "w").write("\n".join(out))
print("wrote", sys.argv[1])
