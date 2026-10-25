#!/usr/bin/env python3
"""Reference values from CoolProp 8 (the reference implementation of
IAPWS-95 and Span et al. 2000) for crates/compressible tests:
single-phase states (T, rho) -> p, u, s, cv, cp, w and saturation states
T -> p, rho_l, rho_v. Usage: python3 build_reference_check.py OUT.csv"""
import sys, CoolProp.CoolProp as CP
rows = []
grids = {
    "Water": ([280, 300, 350, 400, 500, 600, 640, 650, 700, 900, 1200], [0.05, 1, 10, 100, 322, 600, 900, 1000, 1100]),
    "Nitrogen": ([65, 77, 100, 120, 126, 130, 150, 200, 300, 600, 1000], [0.5, 5, 50, 200, 313, 500, 700, 850]),
}
for fl, (ts, rhos) in grids.items():
    for t in ts:
        for rho in rhos:
            try:
                ph = CP.PhaseSI("T", t, "D", rho, fl)
                if "twophase" in ph or ph == "":
                    continue
                vals = [CP.PropsSI(k, "T", t, "D", rho, fl) for k in ("P", "U", "S", "CVMASS", "CPMASS", "A")]
                rows.append(("state", fl, t, rho, *vals))
            except Exception:
                pass
    tt = CP.PropsSI("Ttriple", fl); tc = CP.PropsSI("Tcrit", fl)
    for i in range(1, 20):
        t = tt + (tc - tt) * i / 20.0
        rows.append(("sat", fl, t, 0.0, CP.PropsSI("P", "T", t, "Q", 0, fl), CP.PropsSI("D", "T", t, "Q", 0, fl), CP.PropsSI("D", "T", t, "Q", 1, fl), 0, 0, 0))
with open(sys.argv[1], "w") as f:
    f.write("kind,fluid,t,rho,v1,v2,v3,v4,v5,v6\n")
    for r in rows:
        f.write(",".join(str(x) for x in r) + "\n")
print(len(rows), "rows")
