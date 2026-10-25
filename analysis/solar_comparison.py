"""Plot the fixed GONG transport comparison, without fitting the forecast.

From fluid/: cargo run -p solar --release --example hindcast -- /tmp/solar_forecast.csv
From repo root: python3 analysis/solar_comparison.py /tmp/solar_forecast.csv
Plotting requires numpy and matplotlib; the Rust simulation does not.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def observation(path):
    values = np.loadtxt(path, delimiter=",", comments="#")
    signs = np.where(np.arange(values.shape[1]) % 2 == 0, 1.0, -1.0)
    values -= (values @ signs / values.shape[1])[:, None] * signs
    values -= values.mean()
    return values * 1e4  # tesla -> Gauss


def main():
    root = Path(__file__).resolve().parents[1]
    data = root / "fluid/crates/solar/data"
    initial = observation(data / "gong_20190108.csv")
    target = observation(data / "gong_20190204.csv")
    forecast = np.loadtxt(sys.argv[1], delimiter=",") * 1e4
    if forecast.shape != target.shape or not np.isfinite(forecast).all():
        raise ValueError("forecast must be a finite 60x120 tesla grid")
    fig, axes = plt.subplots(4, 1, figsize=(11, 10), constrained_layout=True)
    titles = ["Initial GONG map: 2019-01-08", "27-day transport: no new emergence",
              "Held-out GONG map: 2019-02-04", "Transport minus held-out map"]
    for ax, values, title in zip(axes, [initial, forecast, target, forecast-target], titles):
        image = ax.imshow(values, origin="lower", extent=[0, 360, -1, 1],
                          aspect="auto", cmap="RdBu_r", vmin=-10, vmax=10,
                          interpolation="nearest")
        ax.set_title(title, loc="left")
        ax.set_ylabel("sin(latitude)")
    axes[-1].set_xlabel("Carrington longitude (degrees)")
    fig.colorbar(image, ax=axes, label="Radial field / difference (Gauss), color clipped at ±10")
    fig.suptitle("Dated synoptic-map comparison; maps span a rotation and use fitted polar fields\n"
                 "Map RMSE: transport 2.936 G; persistence 3.425 G. One case does not establish prediction skill.")
    fig.savefig(root / "analysis/solar_gong_comparison.png", dpi=180)


if __name__ == "__main__":
    main()
