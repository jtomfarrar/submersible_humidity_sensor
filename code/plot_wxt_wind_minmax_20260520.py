"""Plot WXT wind (mean wspd/wdir, wspd_min, wspd_max) from the decoded
20260520 lab-test data produced by uop_decoder.py.

Reads the 1-minute CSV written to
  <WXT_DIR>/20260521_SHS-test_EGH/processed/ASWXT102.csv

Writes PNG to img/Lab20260520/WXT_wind_minmax_20260520.png
"""
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd

REPO    = Path(__file__).resolve().parents[1]
WXTDIR  = Path("/Users/Shared/projects/SHS/20260520/WXT")
CSV     = WXTDIR / "20260521_SHS-test_EGH/processed/ASWXT102.csv"
OUTPNG  = REPO / "img/Lab20260520/WXT_wind_minmax_20260520.png"
OUTPNG.parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    df = pd.read_csv(CSV, index_col="time", parse_dates=True)

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    fig.suptitle(
        f"WXT520 winds — Lab 20260520\n"
        f"{df.index[0]} → {df.index[-1]} UTC   ({len(df)} 1-min records)",
        fontsize=10,
    )

    # ── panel 1: wind speed  mean / min / max ─────────────────────────────
    ax = axes[0]
    ax.fill_between(df.index, df["wspd_min"], df["wspd_max"],
                    color="C0", alpha=0.25, label="wspd_min / max")
    ax.plot(df.index, df["wspd"], "-", lw=1.2, color="C0", label="wspd (vec-avg)")
    ax.set_ylabel("wind speed [m/s]")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, lw=0.4)

    # ── panel 2: wind direction (vector-avg) ──────────────────────────────
    ax = axes[1]
    ax.scatter(df.index, df["wdir"], s=6, color="C1", label="wdir (vec-avg)")
    ax.set_ylabel("wind dir [deg]")
    ax.set_ylim(-5, 365)
    ax.set_yticks([0, 90, 180, 270, 360])
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, lw=0.4)

    # ── panel 3: compass ──────────────────────────────────────────────────
    ax = axes[2]
    ax.scatter(df.index, df["compass"], s=6, color="C4", label="compass (vec-avg)")
    ax.set_ylabel("compass [deg]")
    ax.set_ylim(-5, 365)
    ax.set_yticks([0, 90, 180, 270, 360])
    ax.set_xlabel("time [UTC]")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, lw=0.4)

    for ax in axes:
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTPNG, dpi=140)
    print(f"wrote {OUTPNG}")

    # ── summary ──────────────────────────────────────────────────────────
    print("\n--- wind summary ---")
    print(f"wspd  mean={df['wspd'].mean():.3f}  "
          f"min={df['wspd_min'].min():.3f}  "
          f"max={df['wspd_max'].max():.3f}  m/s")
    print(f"wdir  mean={df['wdir'].mean():.1f}  "
          f"min={df['wdir'].min():.1f}  "
          f"max={df['wdir'].max():.1f}  deg")
    print(f"compass  mean={df['compass'].mean():.1f}  "
          f"min={df['compass'].min():.1f}  "
          f"max={df['compass'].max():.1f}  deg")


if __name__ == "__main__":
    main()
