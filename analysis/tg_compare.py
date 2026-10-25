"""Compare Taylor–Green Re=1600 runs against the 512^3 spectral reference DNS
(HiOCFD workshop, cfd.ku.edu/hiocfd/spectral_Re1600_512.gdiag).

Two dissipation measures are compared, as is standard for this benchmark:
  -dE/dt          (total dissipation, from the energy history)
  2 nu Z          (dissipation computed from enstrophy)
For an exact incompressible solution these agree; their gap measures
under-resolution. Usage: python3 analysis/tg_compare.py analysis/runs/tg_cuda_n*.csv
"""

import csv
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

NU = 1.0 / 1600.0
REF = "fluid/reference/tgv_re1600_spectral512.gdiag"


def load_ref():
    d = np.loadtxt(REF, comments="#")
    return {"t": d[:, 0], "E": d[:, 1], "eps": d[:, 2], "Z": d[:, 3]}


def load_run(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    t = np.array([float(r["t"]) for r in rows])
    E = np.array([float(r["energy"]) for r in rows])
    Z = np.array([float(r["enstrophy"]) for r in rows])
    # central differences for -dE/dt on the uniform output grid
    dEdt = -np.gradient(E, t)
    return {"t": t, "E": E, "Z": Z, "eps_Z": 2 * NU * Z, "eps_E": dEdt}


OUT = "analysis/tg_compare.png"


def main(paths):
    ref = load_ref()
    runs = {}
    print(f"reference: peak -dE/dt = {ref['eps'].max():.6f} at t = {ref['t'][ref['eps'].argmax()]:.2f}")
    print(f"{'run':28s} {'peak eps(2nuZ)':>14s} {'t_peak':>7s} {'max|d eps_E|':>13s} {'max|d eps_Z|':>13s} {'max|dE|/E0':>11s}")
    for p in paths:
        r = load_run(p)
        m = (r["t"] > 0.05) & (r["t"] <= ref["t"].max())
        t = r["t"][m]
        er_E = np.abs(r["eps_E"][m] - np.interp(t, ref["t"], ref["eps"]))
        er_Z = np.abs(r["eps_Z"][m] - np.interp(t, ref["t"], 2 * NU * ref["Z"]))
        er_En = np.abs(r["E"][m] - np.interp(t, ref["t"], ref["E"])) / 0.125
        ip = np.argmax(r["eps_Z"])
        name = p.split("/")[-1]
        print(f"{name:28s} {r['eps_Z'][ip]:14.6f} {r['t'][ip]:7.2f} {er_E.max():13.2e} {er_Z.max():13.2e} {er_En.max():11.2e}")
        runs[name] = r

    BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
    INK, MUTED, GRID, BG = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True, facecolor=BG)
    for ax in (ax1, ax2):
        ax.set_facecolor(BG)
        ax.grid(True, color=GRID, linewidth=0.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    ax1.plot(ref["t"], ref["eps"], color=INK, lw=1.5, ls="--", label="reference 512³ (−dE/dt)")
    colors = [BLUE, ORANGE, AQUA, YELLOW]
    for (name, r), col in zip(runs.items(), colors):
        label = (name.replace("tg_cuda_sph_n", "GPU spherical ").replace("tg_cuda_n", "GPU 2/3-rule ")
                 .replace(".csv", "³"))
        ax1.plot(r["t"], r["eps_Z"], color=col, lw=2, label=f"{label} (2νZ)")
        err = np.abs(r["eps_Z"] - np.interp(r["t"], ref["t"], 2 * NU * ref["Z"]))
        ax2.semilogy(r["t"][1:], np.maximum(err[1:], 1e-12), color=col, lw=2, label=label)
    ax1.set_ylabel("dissipation rate ε")
    ax1.set_title("Taylor–Green vortex, Re = 1600: GPU solver vs. 512³ reference DNS",
                  loc="left", color=INK, fontsize=11)
    ax1.legend(frameon=False, fontsize=9)
    ax2.set_ylabel("|ε − ε_ref|  (enstrophy-based)")
    ax2.set_xlabel("time")
    ax2.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    out = OUT
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0].startswith("--png="):
        OUT = args.pop(0).split("=", 1)[1]
    main(args)
