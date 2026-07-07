"""Plot WXT wind (mean wspd/wdir, wspd_min, wspd_max) from the decoded
20260528 lab-test data produced by uop_decoder.py.

Reads the 1-minute CSV written to
  <WXT_DIR>/20260528_SHS-test_EGH/processed/ASWXT203.csv

Writes PNG to img/Lab20260528/WXT_wind_minmax_20260528.png
"""
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd

REPO    = Path(__file__).resolve().parents[1]
WXTDIR  = Path("/Users/Shared/projects/SHS/20260528/WXT")
CSV     = WXTDIR / "20260528_SHS-test_EGH/processed/ASWXT203.csv"
PKL     = WXTDIR / "20260528_SHS-test_EGH/processed/ASWXT203.pkl"
OUTPNG  = REPO / "img/Lab20260528/WXT_wind_minmax_20260528.png"
OUTPNG.parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    df = pd.read_csv(CSV, index_col="time", parse_dates=True)

    # scalar stats computed from the 11 five-second samples — internally
    # consistent (same population).  wspd_min/wspd_max from the instrument
    # are derived at higher temporal resolution (likely ~1 Hz running averages)
    # and are NOT the min/max of wspd11.
    df_pkl = pd.read_pickle(PKL)
    df["wspd_scalar_mean"] = df_pkl["wspd11"].apply(np.mean)
    df["wspd11_min"]       = df_pkl["wspd11"].apply(np.min)
    df["wspd11_max"]       = df_pkl["wspd11"].apply(np.max)

    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    fig.suptitle(
        f"WXT520 winds — Lab 20260528\n"
        f"{df.index[0]} → {df.index[-1]} UTC   ({len(df)} 1-min records)",
        fontsize=10,
    )

    # ── panel 1: wspd11 min / mean / max  (same 11-sample population) ────
    ax = axes[0]
    ax.fill_between(df.index, df["wspd11_min"], df["wspd11_max"],
                    color="C0", alpha=0.25, label="wspd11 min / max")
    ax.plot(df.index, df["wspd_scalar_mean"], "-", lw=1.2, color="C0",
            label="wspd scalar mean (mean of wspd11)")
    ax.plot(df.index, df["wspd"], "-", lw=1.0, color="C3", alpha=0.7,
            label="wspd vector mean")
    ax.set_ylabel("wind speed [m/s]")
    ax.set_title("Panel 1: from 11×5s samples — min/mean/max are internally consistent", fontsize=8)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, lw=0.4)

    # ── panel 2: instrument wspd_min / wspd_max (high-res, separate quantity)
    ax = axes[1]
    ax.fill_between(df.index, df["wspd_min"], df["wspd_max"],
                    color="C2", alpha=0.25, label="instrument wspd_min / max (high-res)")
    ax.plot(df.index, df["wspd_scalar_mean"], "-", lw=1.0, color="C0", alpha=0.7,
            label="wspd scalar mean (from wspd11, shown for reference)")
    ax.set_ylabel("wind speed [m/s]")
    ax.set_title("Panel 2: instrument-reported min/max (likely ~1 Hz running avg, NOT from wspd11)", fontsize=8)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, lw=0.4)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, lw=0.4)

    # ── panel 3: wind direction (vector-avg) ──────────────────────────────
    ax = axes[2]
    ax.scatter(df.index, df["wdir"], s=4, color="C1", label="wdir (vec-avg)")
    ax.set_ylabel("wind dir [deg]")
    ax.set_ylim(-5, 365)
    ax.set_yticks([0, 90, 180, 270, 360])
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, lw=0.4)

    # ── panel 4: compass ──────────────────────────────────────────────────
    ax = axes[3]
    ax.scatter(df.index, df["compass"], s=4, color="C4", label="compass (vec-avg)")
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
    print(f"wspd11: min={df['wspd11_min'].mean():.3f}  mean={df['wspd_scalar_mean'].mean():.3f}  max={df['wspd11_max'].mean():.3f}  m/s  (mean of per-record stats from 11 samples)")
    print(f"instrument wspd_min mean={df['wspd_min'].mean():.3f}  wspd_max mean={df['wspd_max'].mean():.3f}  m/s  (high-res, separate quantity)")
    print(f"wdir  mean={df['wdir'].mean():.1f}  "
          f"min={df['wdir'].min():.1f}  "
          f"max={df['wdir'].max():.1f}  deg")
    print(f"compass  mean={df['compass'].mean():.1f}  "
          f"min={df['compass'].min():.1f}  "
          f"max={df['compass'].max():.1f}  deg")


if __name__ == "__main__":
    main()
