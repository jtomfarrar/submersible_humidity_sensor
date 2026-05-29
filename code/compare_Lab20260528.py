"""Compare 2026-05-28 SHS lab-test temperatures against the independent
in-tank Seabird GPCTD reference, with WXT520 PTU and attitude diagnostics.

Inputs:
  - data/20260528Lab/RBR/psychrometrics_20260528.csv
  - data/20260528Lab/tank_CT/20260528_tank_CT_log_EGH
  - data/20260528Lab/WXT/processed/ASWXT203.csv

Outputs:
  - data/20260528Lab/compare_SHS_tank_20260528.csv
  - img/Lab20260528/Lab20260528_SHS_vs_tank.png
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vapo_sat_Lab20260528 import DAY, NOTES, note_time


IN_OUT = [
    ("13:11", "13:27", "IN-1"),
    ("13:56", "14:01", "IN-2"),
    ("14:56", "15:01", "IN-3"),
    ("15:42", "15:49", "IN-4"),
    ("16:30", "16:34", "IN-5"),
    ("18:30", "18:50", "IN-6"),
]


def load_tank_ct(path: Path) -> pd.DataFrame:
    """Parse the Seabird GPCTD bench log."""
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
    df = pd.DataFrame(
        {"p_dbar": p, "T_C": T, "C_Sm": C}, index=pd.to_datetime(times)
    )
    df.index.name = "time"
    return df


def main() -> None:
    here = Path(__file__).resolve().parent
    root = here / ".."
    shs_csv = root / "data/20260528Lab/RBR/psychrometrics_20260528.csv"
    ct_log = root / "data/20260528Lab/tank_CT/20260528_tank_CT_log_EGH"
    wxt_csv = root / "data/20260528Lab/WXT/processed/ASWXT203.csv"
    img_dir = root / "img/Lab20260528"
    img_dir.mkdir(parents=True, exist_ok=True)

    shs = pd.read_csv(shs_csv, parse_dates=["time"]).set_index("time")
    print(f"SHS psychrometric samples : {len(shs)}")

    ct = load_tank_ct(ct_log)
    print(
        f"Tank GPCTD samples       : {len(ct)} "
        f"({ct.index.min()} -> {ct.index.max()})"
    )
    ct3 = ct["T_C"].resample("3s").mean()

    wxt = pd.read_csv(wxt_csv, parse_dates=["time"]).set_index("time")
    ptu_dead = bool((wxt[["atmp", "hrh", "bpr"]].abs().sum().sum()) == 0)
    print(
        f"WXT520 records           : {len(wxt)} "
        f"({wxt.index.min()} -> {wxt.index.max()})  "
        f"{'PTU sensor not reporting' if ptu_dead else 'PTU OK'}"
    )

    print("\nWhile SHS is in the tank (independent T reference = GPCTD):")
    print(
        f"{'window':<18s} {'N':>5s} {'T_tank':>8s} "
        f"{'T_dry':>8s} {'T_wet':>8s} "
        f"{'Td-Tt':>8s} {'Tw-Tt':>8s}"
    )
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
        print(
            f"{lbl} {h0}-{h1:<5s} {n:5d} "
            f"{Tt:8.3f} {Td:8.3f} {Tw:8.3f} "
            f"{Td - Tt:+8.3f} {Tw - Tt:+8.3f}"
        )
        rows.append(
            {
                "window": lbl,
                "t0": h0,
                "t1": h1,
                "N": n,
                "T_tank_C": Tt,
                "T_dry_SHS_C": Td,
                "T_wet_SHS_C": Tw,
                "dT_dry_minus_tank": Td - Tt,
                "dT_wet_minus_tank": Tw - Tt,
            }
        )

    summary = pd.DataFrame(rows)
    out_csv = root / "data/20260528Lab/compare_SHS_tank_20260528.csv"
    summary.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    fig, axs = plt.subplots(
        4, 1, figsize=(11, 11), sharex=True, gridspec_kw={"height_ratios": [3, 2, 2, 2]}
    )

    axs[0].plot(shs.index, shs["T_dry_C"], lw=0.7, color="C0", label="SHS T_dry")
    axs[0].plot(shs.index, shs["T_wet_C"], lw=0.7, color="C1", label="SHS T_wet")
    axs[0].plot(ct3.index, ct3.values, lw=1.0, color="C2", label="Tank GPCTD T")
    axs[0].set_ylabel("T (deg C)")
    axs[0].legend(loc="lower right", fontsize=8)
    axs[0].grid(alpha=0.4)
    axs[0].set_title("SHS lab test 2026-05-28 - comparison against in-tank GPCTD")

    diff_dry = np.full(len(shs), np.nan)
    diff_wet = np.full(len(shs), np.nan)
    for h0, h1, _ in IN_OUT:
        t0 = pd.Timestamp(f"{DAY} {h0}")
        t1 = pd.Timestamp(f"{DAY} {h1}")
        sel = (shs.index >= t0) & (shs.index <= t1)
        diff_dry[sel] = shs["T_dry_C"].values[sel] - ct_on_shs.values[sel]
        diff_wet[sel] = shs["T_wet_C"].values[sel] - ct_on_shs.values[sel]
    axs[1].axhline(0, color="0.5", lw=0.5)
    axs[1].plot(shs.index, diff_dry, lw=0.8, color="C0", label="T_dry,SHS - T_tank")
    axs[1].plot(shs.index, diff_wet, lw=0.8, color="C1", label="T_wet,SHS - T_tank")
    axs[1].set_ylabel("Delta T (deg C)")
    axs[1].legend(loc="upper right", fontsize=8)
    axs[1].grid(alpha=0.4)

    axs[2].plot(wxt.index, wxt["atmp"], lw=0.8, color="C3", label="WXT air T")
    axs[2].set_ylabel("WXT T (deg C)", color="C3")
    axs[2].tick_params(axis="y", labelcolor="C3")
    axs[2].grid(alpha=0.4)
    ax2b = axs[2].twinx()
    ax2b.plot(wxt.index, wxt["hrh"], lw=0.8, color="C0", label="WXT RH")
    ax2b.set_ylabel("WXT RH (%)", color="C0")
    ax2b.tick_params(axis="y", labelcolor="C0")
    h1, l1 = axs[2].get_legend_handles_labels()
    h2, l2 = ax2b.get_legend_handles_labels()
    axs[2].legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)

    ax3 = axs[3]
    ax3.plot(wxt.index, wxt["compass"], ".-", ms=2, lw=0.6, color="C4", label="WXT compass")
    ax3.set_ylabel("compass (deg)", color="C4")
    ax3.tick_params(axis="y", labelcolor="C4")
    ax3.grid(alpha=0.4)
    ax3b = ax3.twinx()
    ax3b.plot(wxt.index, wxt["tilt_x"], ".-", ms=2, lw=0.6, color="C5", label="tilt_x")
    ax3b.plot(wxt.index, wxt["tilt_y"], ".-", ms=2, lw=0.6, color="C6", label="tilt_y")
    ax3b.set_ylabel("tilt (deg)", color="0.3")
    h1, l1 = ax3.get_legend_handles_labels()
    h2, l2 = ax3b.get_legend_handles_labels()
    ax3.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)

    color_for = {
        "IN": "tab:blue",
        "OUT": "tab:orange",
        "start/fan-on": "k",
        "OFF": "k",
    }
    for h, kind in NOTES:
        ts = note_time(h)
        c = color_for.get(kind, "gray")
        for ax in axs:
            ax.axvline(ts, color=c, ls="--", lw=0.5, alpha=0.5)
    for h0, h1, _ in IN_OUT:
        t0 = pd.Timestamp(f"{DAY} {h0}")
        t1 = pd.Timestamp(f"{DAY} {h1}")
        for ax in axs:
            ax.axvspan(t0, t1, color="C0", alpha=0.08, zorder=0)

    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axs[-1].set_xlabel("time (UTC)")
    fig.tight_layout()

    out_png = img_dir / "Lab20260528_SHS_vs_tank.png"
    fig.savefig(out_png, dpi=130)
    print(f"saved {out_png}")


if __name__ == "__main__":
    main()
