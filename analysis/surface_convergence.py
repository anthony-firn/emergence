"""Sphere-array drag convergence: point-force walls vs smooth forcing
extension (both GPU), against Sangani & Acrivos K(c) for R = π/3.
Usage: python3 analysis/surface_convergence.py"""

import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAT = re.compile(r"N=\s*(\d+).*err ([+-][\d.]+)%")


def load(path):
    pts = sorted((int(m[1]), abs(float(m[2]))) for m in map(PAT.search, open(path)) if m)
    return [p[0] for p in pts], [p[1] for p in pts]


BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, MUTED, GRID, BG = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(figsize=(7, 4.6), facecolor=BG)
ax.set_facecolor(BG)
ax.grid(True, which="both", color=GRID, linewidth=0.8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for path, col, label in [("analysis/runs/sphere_point_gpu.txt", BLUE, "point-force walls"),
                         ("analysis/runs/sphere_ext_gpu.txt", ORANGE, "smooth forcing extension")]:
    n, e = load(path)
    ax.loglog(n, e, "o-", color=col, lw=2, ms=6, label=label)
    print(label, list(zip(n, e)))
n = [32, 192]
ax.loglog(n, [6.0 * (32 / k) for k in n], ls=":", color=MUTED, lw=1)
ax.text(150, 6.0 * 32 / 150 * 1.3, "1st order", color=MUTED)
ax.set_xlabel("grid size N")
ax.set_ylabel("|K − K_exact| / K_exact  (%)")
ax.set_title("Stokes drag of a periodic sphere array (R = π/3) on the GPU", loc="left", color=INK, fontsize=11)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig("analysis/surface_convergence.png", dpi=150)
print("wrote analysis/surface_convergence.png")
