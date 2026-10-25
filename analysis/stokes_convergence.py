"""Grid convergence of the implicitly penalized sphere against the exact
Stokes drag of a periodic sphere array (Hasimoto 1959; Sangani & Acrivos
1982): K = F/(6πμaU) = 1.2139 at volume fraction c = 1.022e-3.

For each penalization family eta = c_eta·dx²/ν, fit K(N) = K_inf − A·N^(−p)
through the three grids and report the extrapolated K_inf. If the method is
consistent, every family extrapolates to the exact value.

Usage: python3 analysis/stokes_convergence.py analysis/runs/stokes_eta_dx2*.txt
"""

import re
import sys

import numpy as np
from scipy.optimize import brentq

EXACT = 1.2139
PAT = re.compile(r"\[c=([\d.]+)\] N=(\d+) .*K\(momentum\)=([\d.]+)")


def fit(ns, ks):
    """Exact 3-point fit of K = K_inf - A N^-p."""
    (n1, n2, n3), (k1, k2, k3) = ns, ks
    r = (k2 - k1) / (k3 - k2)

    def f(p):
        return (n1**-p - n2**-p) / (n2**-p - n3**-p) - r

    if f(0.02) * f(6.0) > 0:
        return float("nan"), float("nan")  # increments not consistent with a power law
    p = brentq(f, 0.02, 6.0)
    a = (k2 - k1) / (n1**-p - n2**-p)
    return k3 + a * n3**-p, p


def main(paths):
    data = {}
    for path in paths:
        for line in open(path):
            m = PAT.search(line)
            if m:
                c, n, k = float(m[1]), int(m[2]), float(m[3])
                data.setdefault(c, {})[n] = k
    print(f"exact K = {EXACT}")
    print(f"{'eta = c dx^2/nu':>16s} {'K(32)':>8s} {'K(64)':>8s} {'K(96)':>8s} {'order p':>8s} {'K_inf':>8s} {'K_inf err':>9s} {'R1(64,96)':>10s}")
    for c in sorted(data):
        d = data[c]
        ns = sorted(d)
        if len(ns) < 3:
            continue
        ks = [d[n] for n in ns]
        kinf, p = fit(np.array(ns, float), np.array(ks))
        r1 = (ns[2] * ks[2] - ns[1] * ks[1]) / (ns[2] - ns[1])  # first-order Richardson
        print(f"{c:16.2f} {ks[0]:8.4f} {ks[1]:8.4f} {ks[2]:8.4f} {p:8.2f} {kinf:8.4f} {100 * (kinf - EXACT) / EXACT:+8.2f}% {r1:10.4f}")


if __name__ == "__main__":
    main(sys.argv[1:])
