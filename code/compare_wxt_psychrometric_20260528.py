"""Compare WXT520 humidity against SHS psychrometric humidity for 2026-05-28.

The WXT reports air temperature, relative humidity, and pressure.  This script
uses those channels to compute actual vapor pressure and specific humidity,
then compares them against the RBR wet/dry-bulb psychrometric calculation.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vapo_sat_Lab20260528 import DAY, NOTES, es, qs, note_time


REPO = Path(__file__).resolve().parents[1]
WXT_CSV = REPO / "data/20260528Lab/WXT/processed/ASWXT203.csv"
PSY_CSV = REPO / "data/20260528Lab/RBR/psychrometrics_20260528.csv"
OUT_CSV = REPO / "data/20260528Lab/compare_WXT_psychrometric_20260528.csv"
OUT_PNG = REPO / "img/Lab20260528/Lab20260528_WXT_vs_psychrometric_humidity.png"
OUT_PNG.parent.mkdir(parents=True, exist_ok=True)

PSY_PRESSURE_PA = 1013.25 * 100.0

IN_OUT = [
    ("13:11", "13:27", "IN-1"),
    ("13:56", "14:01", "IN-2"),
    ("14:56", "15:01", "IN-3"),
    ("15:42", "15:49", "IN-4"),
    ("16:30", "16:34", "IN-5"),
    ("18:30", "18:50", "IN-6"),
]


def main() -> None:
    wxt = pd.read_csv(WXT_CSV, parse_dates=["time"]).set_index("time")
    psy = pd.read_csv(PSY_CSV, parse_dates=["time"]).set_index("time")

    # Align the 3-second SHS psychrometric series to the 1-minute WXT records.
    psy_on_wxt = psy.reindex(wxt.index, method="nearest", tolerance="5s")
    comp = wxt.join(psy_on_wxt.add_prefix("psy_"))

    wxt_p_pa = comp["bpr"] * 100.0
    comp["wxt_e_hPa"] = es(comp["atmp"].to_numpy(), wxt_p_pa.to_numpy()) * (
        comp["hrh"].to_numpy() / 100.0
    ) / 100.0
    comp["wxt_q_gkg"] = (
        qs(wxt_p_pa.to_numpy(), comp["atmp"].to_numpy())
        * (comp["hrh"].to_numpy() / 100.0)
        * 1000.0
    )
    comp["psy_e_hPa"] = es(comp["psy_T_dry_C"].to_numpy(), PSY_PRESSURE_PA) * (
        comp["psy_RH_pct"].to_numpy() / 100.0
    ) / 100.0
    comp["psy_q_gkg"] = comp["psy_q_kg_per_kg"] * 1000.0
    comp["dRH_psy_minus_wxt_pct"] = comp["psy_RH_pct"] - comp["hrh"]
    comp["de_psy_minus_wxt_hPa"] = comp["psy_e_hPa"] - comp["wxt_e_hPa"]
    comp["dq_psy_minus_wxt_gkg"] = comp["psy_q_gkg"] - comp["wxt_q_gkg"]
    psy_in_air = comp["psy_in_air"].astype("boolean").fillna(False)
    comp["comparison_valid"] = (
        psy_in_air
        & comp[["hrh", "atmp", "bpr", "psy_RH_pct"]].notna().all(axis=1)
        & (comp["wndflag"] == 0)
        & (comp["rhtpflag"] == 0)
    )

    out_cols = [
        "atmp",
        "hrh",
        "bpr",
        "wxt_e_hPa",
        "wxt_q_gkg",
        "psy_T_dry_C",
        "psy_T_wet_C",
        "psy_RH_pct",
        "psy_e_hPa",
        "psy_q_gkg",
        "psy_in_air",
        "comparison_valid",
        "dRH_psy_minus_wxt_pct",
        "de_psy_minus_wxt_hPa",
        "dq_psy_minus_wxt_gkg",
    ]
    comp[out_cols].to_csv(OUT_CSV, index_label="time")
    print(f"wrote {OUT_CSV}")

    valid = comp[comp["comparison_valid"]].copy()
    print_summary("All in-air overlap", valid)

    print("\nAir intervals:")
    air_windows = [
        ("13:08", "13:11", "air-start"),
        ("13:27", "13:56", "air-1"),
        ("14:01", "14:56", "air-2"),
        ("15:01", "15:42", "air-3"),
        ("15:49", "16:30", "air-4"),
        ("16:34", "18:30", "air-5"),
        ("18:50", "2026-05-29 11:08", "air-6 overnight"),
    ]
    for h0, h1, label in air_windows:
        t0 = note_time(h0)
        t1 = note_time(h1)
        sub = valid.loc[(valid.index >= t0) & (valid.index <= t1)]
        print_summary(label, sub, compact=True)

    active_mask = np.zeros(len(valid), dtype=bool)
    for h0, h1, _label in air_windows[1:6]:
        t0 = note_time(h0)
        t1 = note_time(h1)
        active_mask |= (valid.index >= t0) & (valid.index <= t1)
    print_summary("Post-OUT active air", valid.loc[active_mask])

    make_plot(comp, valid)
    print(f"saved {OUT_PNG}")


def print_summary(label: str, df: pd.DataFrame, compact: bool = False) -> None:
    if df.empty:
        print(f"{label:<17s} N=0")
        return
    parts = {
        "N": len(df),
        "WXT_RH": df["hrh"].mean(),
        "PSY_RH": df["psy_RH_pct"].mean(),
        "dRH": df["dRH_psy_minus_wxt_pct"].mean(),
        "WXT_e": df["wxt_e_hPa"].mean(),
        "PSY_e": df["psy_e_hPa"].mean(),
        "de": df["de_psy_minus_wxt_hPa"].mean(),
        "WXT_q": df["wxt_q_gkg"].mean(),
        "PSY_q": df["psy_q_gkg"].mean(),
        "dq": df["dq_psy_minus_wxt_gkg"].mean(),
    }
    if not compact:
        print("\n" + label)
        print(
            "N={N:d}; RH WXT={WXT_RH:.2f}% psych={PSY_RH:.2f}% "
            "diff={dRH:+.2f} pct-pt".format(**parts)
        )
        print(
            "e  WXT={WXT_e:.2f} hPa psych={PSY_e:.2f} hPa "
            "diff={de:+.2f} hPa".format(**parts)
        )
        print(
            "q  WXT={WXT_q:.2f} g/kg psych={PSY_q:.2f} g/kg "
            "diff={dq:+.2f} g/kg".format(**parts)
        )
        return

    print(
        f"{label:<17s} N={parts['N']:4d} "
        f"dRH={parts['dRH']:+6.2f} pct-pt "
        f"de={parts['de']:+6.2f} hPa "
        f"dq={parts['dq']:+6.2f} g/kg"
    )


def make_plot(comp: pd.DataFrame, valid: pd.DataFrame) -> None:
    fig, axs = plt.subplots(4, 1, figsize=(11, 12), sharex=False)

    axs[0].plot(
        comp.index,
        comp["hrh"],
        lw=0.8,
        color="C0",
        label="WXT520 relative humidity (capacitive sensor)",
    )
    axs[0].plot(
        valid.index,
        valid["psy_RH_pct"],
        lw=0.8,
        color="C1",
        label="Psychrometric relative humidity (RBR wet/dry-bulb)",
    )
    axs[0].set_ylabel("Relative humidity (%)")
    axs[0].grid(alpha=0.4)
    axs[0].legend(loc="upper right", fontsize=8)
    axs[0].set_title(
        "2026-05-28 Water vapor: WXT520 capacitive sensor vs RBR psychrometer"
    )

    axs[1].plot(
        comp.index,
        comp["wxt_e_hPa"],
        lw=0.8,
        color="C0",
        label="WXT520 vapor pressure (from T_air, RH, p)",
    )
    axs[1].plot(
        valid.index,
        valid["psy_e_hPa"],
        lw=0.8,
        color="C1",
        label="Psychrometric vapor pressure (from RBR T_dry, T_wet)",
    )
    axs[1].set_ylabel("Actual water vapor pressure, e (hPa)")
    axs[1].grid(alpha=0.4)
    axs[1].legend(loc="upper right", fontsize=8)

    axs[2].axhline(0, color="0.5", lw=0.5)
    axs[2].plot(
        valid.index,
        valid["dRH_psy_minus_wxt_pct"],
        lw=0.8,
        color="C3",
        label="RH difference (psychrometric - WXT520)",
    )
    axs[2].set_ylabel("RH difference (percentage points)")
    axs[2].grid(alpha=0.4)
    ax2 = axs[2].twinx()
    ax2.plot(
        valid.index,
        valid["de_psy_minus_wxt_hPa"],
        lw=0.8,
        color="C4",
        label="Vapor pressure difference (psychrometric - WXT520)",
    )
    ax2.set_ylabel("Vapor pressure difference (hPa)")
    h1, l1 = axs[2].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    axs[2].legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)

    axs[3].plot(
        valid["wxt_e_hPa"],
        valid["psy_e_hPa"],
        ".",
        ms=4,
        alpha=0.55,
        color="C2",
    )
    low = float(np.nanmin([valid["wxt_e_hPa"].min(), valid["psy_e_hPa"].min()]))
    high = float(np.nanmax([valid["wxt_e_hPa"].max(), valid["psy_e_hPa"].max()]))
    axs[3].plot([low, high], [low, high], color="0.4", lw=0.8, label="1:1 line")
    axs[3].set_xlim(low, high)
    axs[3].set_ylim(low, high)
    axs[3].set_xlabel("WXT520 vapor pressure, e (hPa)")
    axs[3].set_ylabel("Psychrometric vapor pressure, e (hPa)")
    axs[3].grid(alpha=0.4)
    axs[3].legend(loc="upper left", fontsize=8)

    for h0, h1, _label in IN_OUT:
        t0 = pd.Timestamp(f"{DAY} {h0}")
        t1 = pd.Timestamp(f"{DAY} {h1}")
        for ax in axs[:3]:
            ax.axvspan(t0, t1, color="0.8", alpha=0.35, zorder=0)

    for h, kind in NOTES:
        ts = note_time(h)
        for ax in axs[:3]:
            ax.axvline(ts, color="0.5", ls="--", lw=0.5, alpha=0.5)

    for ax in axs[:3]:
        ax.set_xlim(comp.index.min(), comp.index.max())
    axs[0].tick_params(labelbottom=False)
    axs[1].tick_params(labelbottom=False)
    axs[2].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axs[2].set_xlabel("time (UTC)")
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130)


if __name__ == "__main__":
    main()
