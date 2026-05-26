"""Compare 2026-05-20 SHS lab-test temperatures against the independent
in-tank Seabird GPCTD reference.

The Vaisala WXT520 file (data/20260520Lab/WXT/ASWXT102.DAT) was checked
separately and its PTU channel decoded to all zeros, so it is not used here.

Inputs:
  - data/20260520Lab/RBR/psychrometrics_20260520.csv (SHS psychrometric output)
  - data/20260520Lab/tank/20260520_CT_tank_log_EGH    (Seabird GPCTD, 1 Hz)

Outputs:
  - data/20260520Lab/compare_SHS_tank_20260520.csv
  - img/Lab20260520/Lab20260520_SHS_vs_tank.png
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from vapo_sat_Lab20260520 import NOTES, DAY


# in-water windows from the notes (start IN -> next OUT)
IN_OUT = [
    ("14:26", "14:33", "IN-1"),
    ("15:07", "15:18", "IN-2"),
    ("15:48", "15:53", "IN-3"),
    ("17:02", "17:08", "IN-4"),
    ("19:18", "19:25", "IN-5"),
]


def load_tank_ct(path: Path) -> pd.DataFrame:
    """Parse the Seabird GPCTD bench log:
        [YYYY-MM-DD HH:MM:SS.fff]  pressure_dbar, temp_C, cond_S/m
    """
    pat = re.compile(
        r"\[(?P<t>[0-9\-]+ [0-9:.]+)\]\s+"
        r"(?P<p>[-0-9.]+),\s+(?P<T>[-0-9.]+),\s+(?P<C>[-0-9.]+)"
    )
    times, p, T, C = [], [], [], []
    with path.open() as f:
        for line in f:
            m = pat.match(line)
            if not m:
                continue
            times.append(m["t"])
            p.append(float(m["p"]))
            T.append(float(m["T"]))
            C.append(float(m["C"]))
    df = pd.DataFrame({"p_dbar": p, "T_C": T, "C_Sm": C},
                      index=pd.to_datetime(times))
    df.index.name = "time"
    return df


def main():
    here = Path(__file__).resolve().parent
    root = here / ".."
    shs_csv = root / "data/20260520Lab/RBR/psychrometrics_20260520.csv"
    ct_log = root / "data/20260520Lab/tank/20260520_CT_tank_log_EGH"
    img_dir = root / "img/Lab20260520"
    img_dir.mkdir(parents=True, exist_ok=True)

    shs = pd.read_csv(shs_csv, parse_dates=["time"]).set_index("time")
    print(f"SHS psychrometric samples : {len(shs)}")

    ct = load_tank_ct(ct_log)
    print(f"Tank GPCTD samples       : {len(ct)} "
          f"({ct.index.min()} → {ct.index.max()})")
    ct3 = ct["T_C"].resample("3s").mean()

    # --- per-in-water-interval comparison -----------------------------
    print("\nWhile SHS is in the tank (independent T reference = GPCTD):")
    print(f"{'window':<18s} {'N':>5s} {'T_tank':>8s} "
          f"{'T_dry':>8s} {'T_wet':>8s} "
          f"{'Td-Tt':>8s} {'Tw-Tt':>8s}")
    rows = []
    ct_on_shs = ct3.reindex(shs.index, method="nearest", tolerance="3s")
    for h0, h1, lbl in IN_OUT:
        t0 = pd.Timestamp(f"{DAY} {h0}")
        t1 = pd.Timestamp(f"{DAY} {h1}")
        s = shs.loc[t0:t1]
        c = ct_on_shs.loc[t0:t1]
        if s.empty or c.dropna().empty:
            continue
        Tt = float(c.mean())
        Td = float(s["T_dry_C"].mean())
        Tw = float(s["T_wet_C"].mean())
        n = int(s["T_dry_C"].notna().sum())
        print(f"{lbl} {h0}-{h1:<5s} {n:5d} "
              f"{Tt:8.3f} {Td:8.3f} {Tw:8.3f} "
              f"{Td-Tt:+8.3f} {Tw-Tt:+8.3f}")
        rows.append({"window": lbl, "t0": h0, "t1": h1, "N": n,
                     "T_tank_C": Tt, "T_dry_SHS_C": Td, "T_wet_SHS_C": Tw,
                     "dT_dry_minus_tank": Td - Tt,
                     "dT_wet_minus_tank": Tw - Tt})

    summary = pd.DataFrame(rows)
    out_csv = root / "data/20260520Lab/compare_SHS_tank_20260520.csv"
    summary.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    # --- plot ---------------------------------------------------------
    fig, axs = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

    axs[0].plot(shs.index, shs["T_dry_C"], lw=0.7, color="C0",
                label="SHS T_dry (RBR)")
    axs[0].plot(shs.index, shs["T_wet_C"], lw=0.7, color="C1",
                label="SHS T_wet (RBR)")
    axs[0].plot(ct3.index, ct3.values, lw=1.0, color="C2",
                label="Tank GPCTD T (3-s avg)")
    axs[0].set_ylabel("T (°C)")
    axs[0].legend(loc="lower right", fontsize=8)
    axs[0].grid(alpha=0.4)
    axs[0].set_title(
        "SHS lab test 2026-05-20 — comparison against in-tank GPCTD"
    )

    # ΔT (SHS − tank) only during in-water windows
    diff_dry = np.full(len(shs), np.nan)
    diff_wet = np.full(len(shs), np.nan)
    for h0, h1, _ in IN_OUT:
        t0 = pd.Timestamp(f"{DAY} {h0}")
        t1 = pd.Timestamp(f"{DAY} {h1}")
        sel = (shs.index >= t0) & (shs.index <= t1)
        diff_dry[sel] = (shs["T_dry_C"].values[sel]
                         - ct_on_shs.values[sel])
        diff_wet[sel] = (shs["T_wet_C"].values[sel]
                         - ct_on_shs.values[sel])
    axs[1].axhline(0, color="0.5", lw=0.5)
    axs[1].plot(shs.index, diff_dry, lw=0.8, color="C0",
                label="T_dry,SHS − T_tank")
    axs[1].plot(shs.index, diff_wet, lw=0.8, color="C1",
                label="T_wet,SHS − T_tank")
    axs[1].set_ylabel("ΔT (°C)")
    axs[1].legend(loc="upper right", fontsize=8)
    axs[1].grid(alpha=0.4)

    color_for = {"IN": "tab:blue", "OUT": "tab:orange",
                 "fan-on/OUT": "tab:red", "start": "k", "OFF": "k"}
    for h, kind in NOTES:
        ts = pd.Timestamp(f"{DAY} {h}")
        c = color_for.get(kind, "gray")
        for ax in axs:
            ax.axvline(ts, color=c, ls="--", lw=0.5, alpha=0.5)
    for h0, h1, _ in IN_OUT:
        t0 = pd.Timestamp(f"{DAY} {h0}")
        t1 = pd.Timestamp(f"{DAY} {h1}")
        for ax in axs:
            ax.axvspan(t0, t1, color="C0", alpha=0.08, zorder=0)

    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[-1].set_xlabel("time (UTC)  2026-05-20")
    fig.tight_layout()

    out_png = img_dir / "Lab20260520_SHS_vs_tank.png"
    fig.savefig(out_png, dpi=130)
    print(f"saved {out_png}")


if __name__ == "__main__":
    main()
