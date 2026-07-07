"""Analyze the 2026-06-29 SHS lab test.

This mirrors the 2026-05-28 lab-test analysis, and keeps the June 29
original-wick test and June 30 plastic-wick follow-up separate where that
matters:

* 2026-06-29: original setup / wick, fan raised to max after first OUT.
* 2026-06-30: plastic wick installed at 12:55 UTC, fan on at first IN.

Inputs:
  - data/20260629Lab/RBR_SHS/233860_20260701_1618.rsk
  - data/20260629Lab/WXT/ASWXT203.DAT or WXT/processed/ASWXT203.csv
  - data/20260629Lab/tank_CT/20260629_tankTC
  - data/20260629Lab/tank_CT/20260630_tankTC.txt

Outputs:
  - data/20260629Lab/RBR_SHS/233860_20260701_1618.nc
  - data/20260629Lab/RBR_SHS/psychrometrics_20260629.csv
  - data/20260629Lab/compare_SHS_tank_20260629.csv
  - data/20260629Lab/compare_WXT_psychrometric_20260629.csv
  - data/20260629Lab/response_time_20260629.csv
  - data/20260629Lab/wet_bulb_settling_20260629.csv
  - img/Lab20260629/*.png
  - docs/Lab20260629/Lab20260629_report.tex
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import curve_fit

from decode_wxt520_sd import decode_wxt520
from rbr_rsk import write_rbr_netcdf


RD_OVER_RV = 287.04 / 461.5
A_FERREL = 0.000662  # Pa^-1 K^-1

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = REPO / "data/20260629Lab"
IMG_DIR = REPO / "img/Lab20260629"
DOC_DIR = REPO / "docs/Lab20260629"

RBR_DIR = DATA_DIR / "RBR_SHS"
RBR_RSK = RBR_DIR / "233860_20260701_1618.rsk"
RBR_NC = RBR_DIR / "233860_20260701_1618.nc"
PSY_CSV = RBR_DIR / "psychrometrics_20260629.csv"

WXT_DAT = DATA_DIR / "WXT/ASWXT203.DAT"
WXT_PROCESSED = DATA_DIR / "WXT/processed"
WXT_CSV = WXT_PROCESSED / "ASWXT203.csv"
WXT_PKL = WXT_PROCESSED / "ASWXT203.pkl"
WXT_META = WXT_PROCESSED / "ASWXT203_meta.json"

CT_LOGS = [
    DATA_DIR / "tank_CT/20260629_tankTC",
    DATA_DIR / "tank_CT/20260630_tankTC.txt",
]

START = pd.Timestamp("2026-06-29 14:50")
END = pd.Timestamp("2026-07-01 16:18")

ACTIVE_WINDOWS = [
    (
        pd.Timestamp("2026-06-29 15:23"),
        pd.Timestamp("2026-06-29 20:24"),
        "2026-06-29 original wick",
    ),
    (
        pd.Timestamp("2026-06-30 12:55"),
        pd.Timestamp("2026-06-30 18:22"),
        "2026-06-30 plastic wick",
    ),
]

IN_WINDOWS = [
    ("2026-06-29 15:28", "2026-06-29 15:33", "J29-IN-1", "2026-06-29 original wick"),
    ("2026-06-29 16:03", "2026-06-29 16:10", "J29-IN-2", "2026-06-29 original wick"),
    ("2026-06-29 16:37", "2026-06-29 16:40", "J29-IN-3", "2026-06-29 original wick"),
    ("2026-06-29 17:18", "2026-06-29 17:23", "J29-IN-4", "2026-06-29 original wick"),
    ("2026-06-29 17:58", "2026-06-29 18:04", "J29-IN-5", "2026-06-29 original wick"),
    ("2026-06-29 18:40", "2026-06-29 18:44", "J29-IN-6", "2026-06-29 original wick"),
    ("2026-06-29 19:15", "2026-06-29 19:20", "J29-IN-7", "2026-06-29 original wick"),
    ("2026-06-29 19:56", "2026-06-29 19:59", "J29-IN-8", "2026-06-29 original wick"),
    ("2026-06-30 13:08", "2026-06-30 13:14", "J30-IN-1", "2026-06-30 plastic wick"),
    ("2026-06-30 13:46", "2026-06-30 13:50", "J30-IN-2", "2026-06-30 plastic wick"),
    ("2026-06-30 14:23", "2026-06-30 14:29", "J30-IN-3", "2026-06-30 plastic wick"),
    ("2026-06-30 14:58", "2026-06-30 15:30", "J30-IN-4", "2026-06-30 plastic wick"),
    ("2026-06-30 16:02", "2026-06-30 16:08", "J30-IN-5", "2026-06-30 plastic wick"),
]

EVENTS = [
    ("2026-06-29 14:50", "start", "RBR SHS start, WXT on"),
    ("2026-06-29 15:23", "fan", "fan level 2"),
    ("2026-06-29 15:27", "tank-log", "tank CT log start"),
    ("2026-06-29 15:28", "IN", "IN"),
    ("2026-06-29 15:33", "OUT", "OUT, fan max"),
    ("2026-06-29 16:03", "IN", "IN"),
    ("2026-06-29 16:10", "OUT", "OUT"),
    ("2026-06-29 16:37", "IN", "IN"),
    ("2026-06-29 16:40", "OUT", "OUT"),
    ("2026-06-29 17:18", "IN", "IN"),
    ("2026-06-29 17:23", "OUT", "OUT"),
    ("2026-06-29 17:58", "IN", "IN"),
    ("2026-06-29 18:04", "OUT", "OUT"),
    ("2026-06-29 18:40", "IN", "IN"),
    ("2026-06-29 18:44", "OUT", "OUT"),
    ("2026-06-29 19:15", "IN", "IN"),
    ("2026-06-29 19:20", "OUT", "OUT"),
    ("2026-06-29 19:56", "IN", "IN"),
    ("2026-06-29 19:59", "OUT", "OUT"),
    ("2026-06-30 12:55", "plastic", "plastic wick, fan off"),
    ("2026-06-30 13:08", "IN", "IN, fan on"),
    ("2026-06-30 13:14", "OUT", "OUT"),
    ("2026-06-30 13:46", "IN", "IN"),
    ("2026-06-30 13:50", "OUT", "OUT"),
    ("2026-06-30 14:23", "IN", "IN"),
    ("2026-06-30 14:29", "OUT", "OUT"),
    ("2026-06-30 14:58", "IN", "IN"),
    ("2026-06-30 15:30", "OUT", "OUT"),
    ("2026-06-30 16:02", "IN", "IN"),
    ("2026-06-30 16:08", "OUT", "OUT"),
]


def es(T_C: np.ndarray, p_Pa: float | np.ndarray) -> np.ndarray:
    """Saturation vapor pressure [Pa] over water (Buck 1981)."""
    p_hPa = np.asarray(p_Pa) * 1e-2
    return 1e2 * 6.1121 * (1.0007 + 3.46e-8 * p_hPa) * np.exp(
        (17.502 * np.asarray(T_C)) / (240.97 + np.asarray(T_C))
    )


def qs(p_Pa: float | np.ndarray, T_C: np.ndarray) -> np.ndarray:
    """Saturation specific humidity [kg/kg]."""
    e = es(T_C, p_Pa)
    return RD_OVER_RV * e / (p_Pa + (RD_OVER_RV - 1.0) * e)


def calc_psy(
    t_dry: np.ndarray, t_wet: np.ndarray, p_Pa: float | np.ndarray
) -> dict[str, np.ndarray]:
    """Return psychrometric depression, RH, dew point, and humidity arrays."""
    depression = t_dry - t_wet
    e_sat_dry = es(t_dry, p_Pa)
    e_sat_wet = es(t_wet, p_Pa)
    e = e_sat_wet - A_FERREL * p_Pa * depression
    rh = 100.0 * e / e_sat_dry
    specific_humidity = qs(p_Pa, t_dry) * rh / 100.0

    dp = t_dry - 10.0
    for _ in range(40):
        e_dp = es(dp, p_Pa)
        e_dp1 = es(dp + 0.1, p_Pa)
        dp = dp + 0.2 * (e - e_dp) / (e_dp1 - e_dp) * 0.1
    return {
        "depression": depression,
        "e": e,
        "rh": rh,
        "dew": dp,
        "specific_humidity": specific_humidity,
    }


def ensure_dirs() -> None:
    for path in (IMG_DIR, DOC_DIR, RBR_DIR, WXT_PROCESSED):
        path.mkdir(parents=True, exist_ok=True)


def ensure_wxt_processed() -> None:
    """Decode ASWXT203.DAT if the processed WXT files are missing."""
    if WXT_CSV.exists() and WXT_PKL.exists() and WXT_META.exists():
        return

    df, meta = decode_wxt520(WXT_DAT)
    df = df.loc[(df.index >= START) & (df.index <= pd.Timestamp("2026-07-01 16:20"))]
    if df.empty:
        raise RuntimeError("No WXT records in the requested analysis interval")
    meta["nrecs_decoded"] = int(len(df))
    meta["time_start"] = df.index.min().isoformat()
    meta["time_end"] = df.index.max().isoformat()

    scalar_cols = [c for c in df.columns if c not in ("wdir11", "wspd11", "compass11")]
    df[scalar_cols].to_csv(WXT_CSV)
    df.to_pickle(WXT_PKL)
    WXT_META.write_text(json.dumps(meta, indent=2))

    ds = xr.Dataset.from_dataframe(df[scalar_cols])
    for col in ("wdir11", "wspd11", "compass11"):
        ds[col] = (("time", "sample"), np.vstack(df[col].to_numpy()))
    ds.to_netcdf(WXT_PROCESSED / "ASWXT203.nc")


def ensure_rbr_nc() -> None:
    if not RBR_NC.exists():
        write_rbr_netcdf(RBR_RSK, RBR_NC)
        print(f"wrote {RBR_NC}")


def events_as_timestamps() -> list[tuple[pd.Timestamp, str, str]]:
    return [(pd.Timestamp(t), kind, label) for t, kind, label in EVENTS]


def state_after(kind: str, current: bool) -> bool:
    if kind == "IN":
        return False
    if kind in {"OUT", "start", "fan", "plastic"}:
        return True
    return current


def in_air_mask(times: pd.DatetimeIndex | np.ndarray) -> np.ndarray:
    t = pd.DatetimeIndex(times)
    mask = np.ones(len(t), dtype=bool)
    current = True
    for ts, kind, _label in events_as_timestamps():
        current = state_after(kind, current)
        mask[t >= ts] = current
    return mask


def phase_for_times(times: pd.DatetimeIndex | np.ndarray) -> np.ndarray:
    t = pd.DatetimeIndex(times)
    phase = np.full(len(t), "", dtype=object)
    for t0, t1, label in ACTIVE_WINDOWS:
        phase[(t >= t0) & (t <= t1)] = label
    return phase


def active_mask(times: pd.DatetimeIndex | np.ndarray) -> np.ndarray:
    return phase_for_times(times) != ""


def load_wxt() -> pd.DataFrame:
    ensure_wxt_processed()
    return pd.read_csv(WXT_CSV, parse_dates=["time"]).set_index("time")


def load_rbr() -> tuple[xr.Dataset, pd.DatetimeIndex, np.ndarray, np.ndarray]:
    ensure_rbr_nc()
    ds = xr.open_dataset(RBR_NC)
    t = pd.DatetimeIndex(pd.to_datetime(ds["time"].values))
    td = ds["temperature"].values.astype(float)
    tw = ds["temperature1"].values.astype(float)
    return ds, t, td, tw


def load_tank_ct_one(path: Path) -> pd.DataFrame:
    pat = re.compile(
        r"\[(?P<t>[0-9\-]+ [0-9:.]+)\]\s+"
        r"(?P<p>[-0-9.]+),\s+(?P<T>[-0-9.]+),\s+(?P<C>[-0-9.]+)"
    )
    rows = []
    with path.open() as f:
        for line in f:
            m = pat.match(line)
            if not m:
                continue
            rows.append(
                {
                    "time": pd.Timestamp(m["t"]),
                    "p_dbar": float(m["p"]),
                    "T_C": float(m["T"]),
                    "C_Sm": float(m["C"]),
                    "source": path.name,
                }
            )
    if not rows:
        raise RuntimeError(f"No tank CT rows parsed from {path}")
    return pd.DataFrame(rows).set_index("time").sort_index()


def load_tank_ct() -> pd.DataFrame:
    frames = [load_tank_ct_one(path) for path in CT_LOGS]
    return pd.concat(frames).sort_index()


def pressure_for_rbr_times(wxt: pd.DataFrame, times: pd.DatetimeIndex, p_default: float) -> np.ndarray:
    wxt_p = wxt["bpr"].dropna()
    if wxt_p.empty:
        return np.full(len(times), p_default)
    return (
        np.interp(
            times.astype(np.int64),
            wxt_p.index.astype(np.int64),
            wxt_p.values,
        )
        * 100.0
    )


def print_interval(
    ts0: pd.Timestamp,
    ts1: pd.Timestamp,
    state: str,
    sel: np.ndarray,
    td: np.ndarray,
    tw: np.ndarray,
    psy: dict[str, np.ndarray],
) -> None:
    n = int(sel.sum())
    if n == 0:
        return
    label = f"{ts0.strftime('%m-%d %H:%M')}-{ts1.strftime('%m-%d %H:%M')}"
    print(
        f"{label:<21s} {state:<6s} {n:5d} "
        f"{np.nanmean(td[sel]):6.2f} "
        f"{np.nanmean(tw[sel]):6.2f} "
        f"{np.nanmean(psy['depression'][sel]):6.2f} "
        f"{np.nanmean(psy['rh'][sel]):6.1f} "
        f"{np.nanmean(psy['dew'][sel]):6.2f} "
        f"{np.nanmean(psy['specific_humidity'][sel]) * 1000:14.2f}"
    )


def psychrometrics() -> pd.DataFrame:
    ds, t_rbr, td, tw = load_rbr()
    wxt = load_wxt()

    p_default = float(ds["default_atmospheric_pressure"].values[0]) * 1e4
    p_pa = pressure_for_rbr_times(wxt, t_rbr, p_default)
    print(
        f"RBR samples              : {len(t_rbr)} "
        f"({t_rbr.min()} -> {t_rbr.max()})"
    )
    print(
        f"using WXT pressure       : "
        f"min={np.nanmin(p_pa)/100:.2f} mean={np.nanmean(p_pa)/100:.2f} "
        f"max={np.nanmax(p_pa)/100:.2f} hPa"
    )

    psy = calc_psy(td, tw, p_pa)
    air = in_air_mask(t_rbr)
    phase = phase_for_times(t_rbr)

    print("\nMean psychrometric values for note-defined intervals:\n")
    print(
        f"{'interval':<21s} {'state':<6s} {'N':>5s} "
        f"{'Td':>6s} {'Tw':>6s} {'dep':>6s} "
        f"{'RH%':>6s} {'Tdew':>6s} {'spec_hum g/kg':>14s}"
    )

    event_ts = events_as_timestamps()
    interval_start = pd.Timestamp(t_rbr[0])
    state = "air"
    t_last = pd.Timestamp(t_rbr[-1])
    for ts, kind, _label in event_ts:
        if ts <= interval_start:
            state = "air" if state_after(kind, state == "air") else "water"
            continue
        if ts > t_last:
            break
        sel = (t_rbr >= interval_start) & (t_rbr < ts)
        print_interval(interval_start, ts, state, sel, td, tw, psy)
        state = "air" if state_after(kind, state == "air") else "water"
        interval_start = ts
    sel = t_rbr >= interval_start
    print_interval(interval_start, t_last, state, sel, td, tw, psy)

    df = pd.DataFrame(
        {
            "time": t_rbr,
            "T_dry_C": td,
            "T_wet_C": tw,
            "depression_C": psy["depression"],
            "RH_pct": psy["rh"],
            "dew_point_C": psy["dew"],
            "specific_humidity_kg_per_kg": psy["specific_humidity"],
            "pressure_used_hPa": p_pa / 100.0,
            "in_air": air,
            "phase": phase,
        }
    )
    df.to_csv(PSY_CSV, index=False)
    print(f"\nwrote {PSY_CSV}")

    plot_psychrometrics(df, wxt, p_default / 100.0)
    return df


def plot_psychrometrics(psy_df: pd.DataFrame, wxt: pd.DataFrame, p_default_hpa: float) -> None:
    df = psy_df.set_index("time")
    fig, axs = plt.subplots(5, 1, figsize=(12, 13), sharex=True)

    axs[0].plot(df.index, df["T_dry_C"], lw=0.7, label="T_dry (RBR temperature)")
    axs[0].plot(df.index, df["T_wet_C"], lw=0.7, label="T_wet (RBR temperature1)")
    axs[0].plot(wxt.index, wxt["atmp"], lw=0.8, color="C3", label="T_air (WXT520 atmp)")
    axs[0].set_ylabel("T (deg C)")
    axs[0].legend(loc="best", fontsize=8)
    axs[0].grid(alpha=0.4)
    axs[0].set_title("Lab20260629 SHS psychrometric humidity")

    axs[1].plot(df.index, df["depression_C"], lw=0.7, color="C3")
    axs[1].set_ylabel("Td - Tw (deg C)")
    axs[1].grid(alpha=0.4)

    rh_air = df["RH_pct"].where(df["in_air"])
    rh_sub = df["RH_pct"].where(~df["in_air"])
    axs[2].plot(df.index, rh_air, lw=0.7, color="C0", label="RH (in air)")
    axs[2].plot(df.index, rh_sub, lw=0.7, color="0.7", label="RH (submerged, not physical)")
    axs[2].plot(wxt.index, wxt["hrh"], lw=0.8, color="C3", alpha=0.7, label="WXT520 RH")
    axs[2].set_ylabel("RH (%)")
    axs[2].set_ylim(0, 115)
    axs[2].legend(loc="best", fontsize=8)
    axs[2].grid(alpha=0.4)

    axs[3].plot(
        df.index,
        df["specific_humidity_kg_per_kg"].where(df["in_air"]) * 1000.0,
        lw=0.7,
        color="C2",
    )
    axs[3].set_ylabel("specific humidity (g/kg)")
    axs[3].grid(alpha=0.4)

    axs[4].plot(df.index, df["pressure_used_hPa"], lw=0.7, color="C0", label="pressure used")
    axs[4].plot(wxt.index, wxt["bpr"], lw=0.8, color="C2", alpha=0.6, label="WXT bpr")
    axs[4].axhline(p_default_hpa, color="C3", lw=1.0, ls="--", label="RBR default pressure")
    axs[4].set_ylabel("pressure (hPa)")
    axs[4].legend(loc="best", fontsize=8)
    axs[4].grid(alpha=0.4)

    mark_events(axs, annotate_axis=axs[0])
    for t0, t1, _label in ACTIVE_WINDOWS:
        for ax in axs:
            ax.axvspan(t0, t1, color="0.9", alpha=0.35, zorder=0)

    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axs[-1].set_xlabel("time (UTC)")
    fig.tight_layout()
    out = IMG_DIR / "Lab20260629_SHS_humidity.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"saved {out}")


def mark_events(axs: list[plt.Axes] | np.ndarray, annotate_axis: plt.Axes | None = None) -> None:
    color_for = {
        "IN": "tab:blue",
        "OUT": "tab:orange",
        "start": "k",
        "fan": "0.2",
        "plastic": "tab:green",
        "tank-log": "0.5",
    }
    for ts, kind, label in events_as_timestamps():
        color = color_for.get(kind, "0.5")
        for ax in axs:
            ax.axvline(ts, color=color, ls="--", lw=0.5, alpha=0.45)
        if annotate_axis is not None and kind in {"IN", "OUT", "plastic"}:
            annotate_axis.annotate(
                label,
                xy=(ts, annotate_axis.get_ylim()[1]),
                xytext=(0, -2),
                textcoords="offset points",
                rotation=90,
                ha="right",
                va="top",
                fontsize=6,
                color=color,
            )


def shade_in_water(axs: list[plt.Axes] | np.ndarray, alpha: float = 0.28) -> None:
    for t0_s, t1_s, _label, _phase in IN_WINDOWS:
        t0 = pd.Timestamp(t0_s)
        t1 = pd.Timestamp(t1_s)
        for ax in np.ravel(axs):
            ax.axvspan(t0, t1, color="0.82", alpha=alpha, zorder=0)


def mark_configuration_change(axs: list[plt.Axes] | np.ndarray) -> None:
    ts = pd.Timestamp("2026-06-30 12:55")
    for ax in np.ravel(axs):
        ax.axvline(ts, color="tab:green", ls="--", lw=0.8, alpha=0.8)


def compare_shs_tank(psy_df: pd.DataFrame) -> pd.DataFrame:
    shs = psy_df.set_index("time")
    ct = load_tank_ct()
    ct3 = ct["T_C"].resample("3s").mean()
    ct_on_shs = ct3.reindex(shs.index, method="nearest", tolerance="3s")
    wxt = load_wxt()

    print(
        f"\nTank CT samples          : {len(ct)} "
        f"({ct.index.min()} -> {ct.index.max()})"
    )
    print("\nWhile SHS is in the tank (independent T reference = tank CT):")
    print(
        f"{'window':<24s} {'N':>5s} {'T_tank':>8s} "
        f"{'T_dry':>8s} {'T_wet':>8s} {'Td-Tt':>8s} {'Tw-Tt':>8s}"
    )

    rows = []
    for t0_s, t1_s, label, phase in IN_WINDOWS:
        t0 = pd.Timestamp(t0_s)
        t1 = pd.Timestamp(t1_s)
        s = shs.loc[t0:t1]
        c = ct_on_shs.loc[t0:t1]
        if s.empty or c.dropna().empty:
            continue
        tank_t = float(c.mean())
        dry_t = float(s["T_dry_C"].mean())
        wet_t = float(s["T_wet_C"].mean())
        n = int(s["T_dry_C"].notna().sum())
        print(
            f"{label} {t0:%m-%d %H:%M}-{t1:%H:%M} {n:5d} "
            f"{tank_t:8.3f} {dry_t:8.3f} {wet_t:8.3f} "
            f"{dry_t - tank_t:+8.3f} {wet_t - tank_t:+8.3f}"
        )
        rows.append(
            {
                "window": label,
                "phase": phase,
                "t0": t0,
                "t1": t1,
                "N": n,
                "T_tank_C": tank_t,
                "T_dry_SHS_C": dry_t,
                "T_wet_SHS_C": wet_t,
                "dT_dry_minus_tank": dry_t - tank_t,
                "dT_wet_minus_tank": wet_t - tank_t,
            }
        )

    summary = pd.DataFrame(rows)
    out_csv = DATA_DIR / "compare_SHS_tank_20260629.csv"
    summary.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    fig, axs = plt.subplots(
        4, 1, figsize=(12, 11), sharex=True, gridspec_kw={"height_ratios": [3, 2, 2, 2]}
    )
    axs[0].plot(shs.index, shs["T_dry_C"], lw=0.7, color="C0", label="SHS T_dry")
    axs[0].plot(shs.index, shs["T_wet_C"], lw=0.7, color="C1", label="SHS T_wet")
    axs[0].plot(ct3.index, ct3.values, lw=0.7, color="C2", label="Tank CT T")
    axs[0].set_ylabel("T (deg C)")
    axs[0].legend(loc="best", fontsize=8)
    axs[0].grid(alpha=0.4)
    axs[0].set_title("Lab20260629 SHS temperatures vs in-tank CT")

    diff_dry = pd.Series(np.nan, index=shs.index)
    diff_wet = pd.Series(np.nan, index=shs.index)
    for t0_s, t1_s, _label, _phase in IN_WINDOWS:
        t0 = pd.Timestamp(t0_s)
        t1 = pd.Timestamp(t1_s)
        sel = (shs.index >= t0) & (shs.index <= t1)
        diff_dry.loc[sel] = shs.loc[sel, "T_dry_C"] - ct_on_shs.loc[sel]
        diff_wet.loc[sel] = shs.loc[sel, "T_wet_C"] - ct_on_shs.loc[sel]
    axs[1].axhline(0, color="0.5", lw=0.5)
    axs[1].plot(diff_dry.index, diff_dry.values, lw=0.8, color="C0", label="T_dry,SHS - T_tank")
    axs[1].plot(diff_wet.index, diff_wet.values, lw=0.8, color="C1", label="T_wet,SHS - T_tank")
    axs[1].set_ylabel("delta T (deg C)")
    axs[1].legend(loc="best", fontsize=8)
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
    axs[2].legend(h1 + h2, l1 + l2, loc="best", fontsize=8)

    axs[3].plot(wxt.index, wxt["compass"], ".-", ms=2, lw=0.5, color="C4", label="WXT compass")
    axs[3].set_ylabel("compass (deg)", color="C4")
    axs[3].tick_params(axis="y", labelcolor="C4")
    axs[3].grid(alpha=0.4)
    ax3b = axs[3].twinx()
    ax3b.plot(wxt.index, wxt["tilt_x"], ".-", ms=2, lw=0.5, color="C5", label="tilt_x")
    ax3b.plot(wxt.index, wxt["tilt_y"], ".-", ms=2, lw=0.5, color="C6", label="tilt_y")
    ax3b.set_ylabel("tilt (deg)", color="0.3")
    h1, l1 = axs[3].get_legend_handles_labels()
    h2, l2 = ax3b.get_legend_handles_labels()
    axs[3].legend(h1 + h2, l1 + l2, loc="best", fontsize=8)

    for t0_s, t1_s, _label, _phase in IN_WINDOWS:
        t0 = pd.Timestamp(t0_s)
        t1 = pd.Timestamp(t1_s)
        for ax in axs:
            ax.axvspan(t0, t1, color="C0", alpha=0.08, zorder=0)
    mark_events(axs)

    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axs[-1].set_xlabel("time (UTC)")
    fig.tight_layout()
    out_png = IMG_DIR / "Lab20260629_SHS_vs_tank.png"
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f"saved {out_png}")

    make_tank_report_plot(shs, ct3, diff_dry, diff_wet)
    return summary


def make_tank_report_plot(
    shs: pd.DataFrame,
    ct3: pd.Series,
    diff_dry: pd.Series,
    diff_wet: pd.Series,
) -> None:
    fig, axs = plt.subplots(2, 1, figsize=(11, 6.8), sharex=True)

    axs[0].plot(shs.index, shs["T_dry_C"], lw=0.8, color="C0", label="SHS T_dry")
    axs[0].plot(shs.index, shs["T_wet_C"], lw=0.8, color="C1", label="SHS T_wet")
    axs[0].plot(ct3.index, ct3.values, lw=0.9, color="C2", label="Tank CT T")
    axs[0].set_ylabel("T (deg C)")
    axs[0].legend(loc="best", fontsize=8)
    axs[0].grid(alpha=0.35)
    axs[0].set_title("SHS temperatures against in-tank CT reference")

    axs[1].axhline(0, color="0.5", lw=0.6)
    axs[1].plot(diff_dry.index, diff_dry.values, lw=0.8, color="C0", label="T_dry - T_tank")
    axs[1].plot(diff_wet.index, diff_wet.values, lw=0.8, color="C1", label="T_wet - T_tank")
    axs[1].set_ylabel("delta T (deg C)")
    axs[1].legend(loc="best", fontsize=8)
    axs[1].grid(alpha=0.35)

    shade_in_water(axs, alpha=0.22)
    mark_configuration_change(axs)
    axs[-1].set_xlim(pd.Timestamp("2026-06-29 15:20"), pd.Timestamp("2026-06-30 16:45"))
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axs[-1].set_xlabel("time (UTC)")
    fig.tight_layout()
    out_png = IMG_DIR / "Lab20260629_SHS_vs_tank_report.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"saved {out_png}")


def compare_wxt_psychrometric(psy_df: pd.DataFrame) -> pd.DataFrame:
    wxt = load_wxt()
    psy = psy_df.set_index("time")
    psy_on_wxt = psy.reindex(wxt.index, method="nearest", tolerance="5s")
    comp = wxt.join(psy_on_wxt.add_prefix("psy_"))

    wxt_p_pa = comp["bpr"] * 100.0
    comp["wxt_e_hPa"] = es(comp["atmp"].to_numpy(), wxt_p_pa.to_numpy()) * (
        comp["hrh"].to_numpy() / 100.0
    ) / 100.0
    comp["wxt_specific_humidity_gkg"] = (
        qs(wxt_p_pa.to_numpy(), comp["atmp"].to_numpy())
        * (comp["hrh"].to_numpy() / 100.0)
        * 1000.0
    )
    comp["psy_e_hPa"] = es(
        comp["psy_T_dry_C"].to_numpy(),
        comp["psy_pressure_used_hPa"].to_numpy() * 100.0,
    ) * (comp["psy_RH_pct"].to_numpy() / 100.0) / 100.0
    comp["psy_specific_humidity_gkg"] = comp["psy_specific_humidity_kg_per_kg"] * 1000.0
    comp["dRH_psy_minus_wxt_pct"] = comp["psy_RH_pct"] - comp["hrh"]
    comp["de_psy_minus_wxt_hPa"] = comp["psy_e_hPa"] - comp["wxt_e_hPa"]
    comp["dspecific_humidity_psy_minus_wxt_gkg"] = (
        comp["psy_specific_humidity_gkg"] - comp["wxt_specific_humidity_gkg"]
    )
    psy_in_air = comp["psy_in_air"].astype("boolean").fillna(False)
    comp["phase"] = phase_for_times(comp.index)
    comp["active_phase"] = comp["phase"] != ""
    comp["comparison_valid"] = (
        psy_in_air
        & comp[["hrh", "atmp", "bpr", "psy_RH_pct"]].notna().all(axis=1)
        & (comp["wndflag"] == 0)
        & (comp["rhtpflag"] == 0)
    )
    comp["active_comparison_valid"] = comp["comparison_valid"] & comp["active_phase"]

    out_cols = [
        "phase",
        "active_phase",
        "atmp",
        "hrh",
        "bpr",
        "wxt_e_hPa",
        "wxt_specific_humidity_gkg",
        "psy_T_dry_C",
        "psy_T_wet_C",
        "psy_RH_pct",
        "psy_e_hPa",
        "psy_specific_humidity_gkg",
        "psy_in_air",
        "comparison_valid",
        "active_comparison_valid",
        "dRH_psy_minus_wxt_pct",
        "de_psy_minus_wxt_hPa",
        "dspecific_humidity_psy_minus_wxt_gkg",
    ]
    out_csv = DATA_DIR / "compare_WXT_psychrometric_20260629.csv"
    comp[out_cols].to_csv(out_csv, index_label="time")
    print(f"\nwrote {out_csv}")

    valid = comp[comp["comparison_valid"]].copy()
    active_valid = comp[comp["active_comparison_valid"]].copy()
    print_summary("All in-air overlap", valid)
    print_summary("Active in-air overlap", active_valid)
    for _t0, _t1, phase in ACTIVE_WINDOWS:
        print_summary(phase, active_valid[active_valid["phase"] == phase], compact=True)

    print("\nAir intervals:")
    for t0, t1, label, phase in air_windows():
        sub = valid.loc[(valid.index >= t0) & (valid.index <= t1)]
        print_summary(f"{label}", sub, compact=True)

    make_wxt_psy_plot(comp, valid)
    make_wxt_psy_report_plot(comp, valid)
    return comp


def print_summary(label: str, df: pd.DataFrame, compact: bool = False) -> None:
    if df.empty:
        print(f"{label:<28s} N=0")
        return
    parts = {
        "N": len(df),
        "WXT_RH": df["hrh"].mean(),
        "PSY_RH": df["psy_RH_pct"].mean(),
        "dRH": df["dRH_psy_minus_wxt_pct"].mean(),
        "WXT_e": df["wxt_e_hPa"].mean(),
        "PSY_e": df["psy_e_hPa"].mean(),
        "de": df["de_psy_minus_wxt_hPa"].mean(),
        "WXT_spec_hum": df["wxt_specific_humidity_gkg"].mean(),
        "PSY_spec_hum": df["psy_specific_humidity_gkg"].mean(),
        "d_spec_hum": df["dspecific_humidity_psy_minus_wxt_gkg"].mean(),
    }
    if compact:
        print(
            f"{label:<28s} N={parts['N']:4d} "
            f"dRH={parts['dRH']:+7.2f} pct-pt "
            f"de={parts['de']:+7.2f} hPa "
            f"d_spec_hum={parts['d_spec_hum']:+7.2f} g/kg"
        )
        return

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
        "specific humidity WXT={WXT_spec_hum:.2f} g/kg "
        "psych={PSY_spec_hum:.2f} g/kg "
        "diff={d_spec_hum:+.2f} g/kg".format(**parts)
    )


def make_wxt_psy_plot(comp: pd.DataFrame, valid: pd.DataFrame) -> None:
    fig, axs = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

    axs[0].plot(comp.index, comp["hrh"], lw=0.8, color="C0", label="WXT520 RH")
    axs[0].plot(valid.index, valid["psy_RH_pct"], lw=0.8, color="C1", label="Psychrometric RH")
    axs[0].set_ylabel("RH (%)")
    axs[0].grid(alpha=0.4)
    axs[0].legend(loc="best", fontsize=8)
    axs[0].set_title("Lab20260629 water vapor: WXT520 vs RBR psychrometer")

    axs[1].plot(comp.index, comp["wxt_e_hPa"], lw=0.8, color="C0", label="WXT520 e")
    axs[1].plot(valid.index, valid["psy_e_hPa"], lw=0.8, color="C1", label="Psychrometric e")
    axs[1].set_ylabel("e (hPa)")
    axs[1].grid(alpha=0.4)
    axs[1].legend(loc="best", fontsize=8)

    axs[2].plot(
        comp.index,
        comp["wxt_specific_humidity_gkg"],
        lw=0.8,
        color="C0",
        label="WXT specific humidity",
    )
    axs[2].plot(
        valid.index,
        valid["psy_specific_humidity_gkg"],
        lw=0.8,
        color="C1",
        label="Psychrometric specific humidity",
    )
    axs[2].set_ylabel("specific humidity (g/kg)")
    axs[2].grid(alpha=0.4)
    axs[2].legend(loc="best", fontsize=8)

    axs[3].axhline(0, color="0.5", lw=0.5)
    axs[3].plot(
        valid.index,
        valid["dRH_psy_minus_wxt_pct"],
        lw=0.8,
        color="C3",
        label="RH diff (psych - WXT)",
    )
    axs[3].set_ylabel("RH diff (pct-pt)")
    axs[3].grid(alpha=0.4)
    ax2 = axs[3].twinx()
    ax2.plot(
        valid.index,
        valid["de_psy_minus_wxt_hPa"],
        lw=0.8,
        color="C4",
        label="e diff (psych - WXT)",
    )
    ax2.set_ylabel("e diff (hPa)")
    h1, l1 = axs[3].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    axs[3].legend(h1 + h2, l1 + l2, loc="best", fontsize=8)

    for t0_s, t1_s, _label, _phase in IN_WINDOWS:
        t0 = pd.Timestamp(t0_s)
        t1 = pd.Timestamp(t1_s)
        for ax in axs:
            ax.axvspan(t0, t1, color="0.8", alpha=0.3, zorder=0)
    mark_events(axs)

    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axs[-1].set_xlabel("time (UTC)")
    fig.tight_layout()
    out_png = IMG_DIR / "Lab20260629_WXT_vs_psychrometric_humidity.png"
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f"saved {out_png}")


def make_wxt_psy_report_plot(comp: pd.DataFrame, valid: pd.DataFrame) -> None:
    fig, axs = plt.subplots(2, 1, figsize=(11, 6.8), sharex=True)

    axs[0].plot(comp.index, comp["hrh"], lw=0.8, color="C0", label="WXT520 RH")
    axs[0].plot(valid.index, valid["psy_RH_pct"], lw=0.8, color="C1", label="Psychrometric RH")
    axs[0].set_ylabel("RH (%)")
    axs[0].legend(loc="best", fontsize=8)
    axs[0].grid(alpha=0.35)
    axs[0].set_title("WXT520 humidity vs RBR wet/dry-bulb psychrometer")

    axs[1].axhline(0, color="0.5", lw=0.6)
    axs[1].plot(
        valid.index,
        valid["dRH_psy_minus_wxt_pct"],
        lw=0.8,
        color="C3",
        label="RH diff (psych - WXT)",
    )
    axs[1].set_ylabel("RH diff (pct-pt)")
    axs[1].grid(alpha=0.35)
    ax2 = axs[1].twinx()
    ax2.plot(
        valid.index,
        valid["de_psy_minus_wxt_hPa"],
        lw=0.8,
        color="C4",
        label="e diff (psych - WXT)",
    )
    ax2.set_ylabel("e diff (hPa)")
    h1, l1 = axs[1].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    axs[1].legend(h1 + h2, l1 + l2, loc="best", fontsize=8)

    shade_in_water(axs, alpha=0.22)
    mark_configuration_change(axs)
    axs[-1].set_xlim(pd.Timestamp("2026-06-29 15:00"), END)
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axs[-1].set_xlabel("time (UTC)")
    fig.tight_layout()
    out_png = IMG_DIR / "Lab20260629_WXT_vs_psychrometric_humidity_report.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"saved {out_png}")


def air_windows() -> list[tuple[pd.Timestamp, pd.Timestamp, str, str]]:
    windows = []
    for i, (ts, kind, _label) in enumerate(events_as_timestamps()):
        if kind != "OUT":
            continue
        ts_end = None
        for ts2, kind2, _label2 in events_as_timestamps()[i + 1 :]:
            if kind2 == "IN":
                ts_end = ts2
                break
        max_end = min(ts + pd.Timedelta(minutes=25), END)
        if ts_end is None or ts_end > max_end:
            ts_end = max_end
        phase = phase_for_times([ts])[0]
        label = f"{ts:%m-%d %H:%M}-{ts_end:%H:%M}"
        windows.append((ts, ts_end, label, phase))
    # Include pre-dunk fan periods, mostly for table completeness.
    windows.insert(
        0,
        (
            pd.Timestamp("2026-06-29 15:23"),
            pd.Timestamp("2026-06-29 15:28"),
            "06-29 15:23-15:28",
            "2026-06-29 original wick",
        ),
    )
    windows.insert(
        9,
        (
            pd.Timestamp("2026-06-30 12:55"),
            pd.Timestamp("2026-06-30 13:08"),
            "06-30 12:55-13:08",
            "2026-06-30 plastic wick",
        ),
    )
    return windows


def _fit_tau(t_sec: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray]:
    """Fit y = A + B exp(-t/tau). Return tau_s and fitted y."""
    ok = np.isfinite(t_sec) & np.isfinite(y)
    t_sec = t_sec[ok]
    y = y[ok]
    if len(t_sec) < 8 or np.nanmax(y) - np.nanmin(y) < 1e-6:
        return float("nan"), np.full(len(ok), np.nan)

    def _model(t, a, b, tau):
        return a + b * np.exp(-t / tau)

    a0 = float(y[-1])
    b0 = float(y[0] - y[-1])
    tau0 = max(float(t_sec[-1]) / 3.0, 10.0)
    try:
        popt, _ = curve_fit(
            _model,
            t_sec,
            y,
            p0=[a0, b0, tau0],
            bounds=([-np.inf, -np.inf, 1.0], [np.inf, np.inf, 7200.0]),
            maxfev=10000,
        )
        y_fit_valid = _model(t_sec, *popt)
        y_fit = np.full(len(ok), np.nan)
        y_fit[np.where(ok)[0]] = y_fit_valid
        return float(popt[2]), y_fit
    except (RuntimeError, ValueError):
        return float("nan"), np.full(len(ok), np.nan)


def response_times(psy_df: pd.DataFrame) -> pd.DataFrame:
    df = psy_df.set_index("time")
    out_events = []
    ev = events_as_timestamps()
    for i, (ts, kind, _label) in enumerate(ev):
        if kind != "OUT":
            continue
        next_in = None
        for ts2, kind2, _label2 in ev[i + 1 :]:
            if kind2 == "IN":
                next_in = ts2
                break
        fit_end = ts + pd.Timedelta(minutes=25)
        if next_in is not None and next_in < fit_end:
            fit_end = next_in
        out_events.append((ts, fit_end, phase_for_times([ts])[0]))

    var_specs = [
        ("T_dry", "T_dry_C", "T dry-bulb (deg C)"),
        ("T_wet", "T_wet_C", "T wet-bulb (deg C)"),
        ("depression", "depression_C", "Depression Td - Tw (deg C)"),
        ("RH_psychro", "RH_pct", "Psychrometric RH (%)"),
    ]

    all_fits: dict[str, list] = {name: [] for name, _col, _ylabel in var_specs}
    rows = []

    print("\n--- response-time estimates (all OUT events, <=25-min windows) ---")
    print(
        f"  {'event':<18s} {'phase':<14s} {'T_dry':>11s} {'T_wet':>11s} "
        f"{'depression':>11s} {'RH_psychro':>11s}"
    )
    for t_out, fit_end, phase in out_events:
        seg = df.loc[t_out:fit_end]
        if len(seg) < 8:
            continue
        t_sec = (seg.index - t_out).total_seconds().astype(float).to_numpy()
        row = {
            "event": f"OUT {t_out:%m-%d %H:%M}",
            "t_out": t_out,
            "t_fit_end": fit_end,
            "phase": phase,
            "fit_window_min": (fit_end - t_out).total_seconds() / 60.0,
        }
        taus = []
        for name, col, _ylabel in var_specs:
            tau, y_fit = _fit_tau(t_sec, seg[col].to_numpy(dtype=float))
            all_fits[name].append((t_out, fit_end, seg.index, seg[col].to_numpy(dtype=float), y_fit, tau))
            row[f"tau_{name}_s"] = tau
            taus.append(tau)
        rows.append(row)
        tau_str = "  ".join("--" if not np.isfinite(x) else f"{x:6.0f} s" for x in taus)
        print(f"  OUT {t_out:%m-%d %H:%M} {phase[:14]:<14s} {tau_str}")

    summary = pd.DataFrame(rows)
    out_csv = DATA_DIR / "response_time_20260629.csv"
    summary.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    fig, axs = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, (name, _col, ylabel) in zip(axs, var_specs):
        for i, (t_out, fit_end, idx, y_data, y_fit, tau) in enumerate(all_fits[name]):
            color = colors[i % len(colors)]
            t_min = (idx - t_out).total_seconds() / 60.0
            ax.plot(t_min, y_data, ".", ms=2.5, color=color, alpha=0.55)
            label = f"OUT {t_out:%m-%d %H:%M}, tau={tau:.0f}s" if np.isfinite(tau) else f"OUT {t_out:%m-%d %H:%M}"
            ax.plot(t_min, y_fit, "-", lw=1.2, color=color, label=label)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(alpha=0.4)
        ax.legend(loc="best", fontsize=6, ncol=2)
    axs[0].set_title("RBRduo3 response after each OUT event - Lab20260629")
    axs[-1].set_xlabel("time since OUT (min)")
    fig.tight_layout()
    out_png = IMG_DIR / "Lab20260629_response_time.png"
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    print(f"saved {out_png}")

    make_response_report_plot(all_fits)
    return summary


def make_response_report_plot(all_fits: dict[str, list]) -> None:
    specs = [
        ("T_wet", "T wet-bulb (deg C)"),
        ("RH_psychro", "Psychrometric RH (%)"),
    ]
    fig, axs = plt.subplots(2, 1, figsize=(11, 7.0), sharex=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, (name, ylabel) in zip(axs, specs):
        for i, (t_out, _fit_end, idx, y_data, y_fit, tau) in enumerate(all_fits[name]):
            color = colors[i % len(colors)]
            t_min = (idx - t_out).total_seconds() / 60.0
            ax.plot(t_min, y_data, ".", ms=2.2, color=color, alpha=0.45)
            if np.isfinite(tau) and tau < 0.98 * 7200.0:
                label = f"OUT {t_out:%m-%d %H:%M}, tau={tau:.0f}s"
            elif np.isfinite(tau):
                label = f"OUT {t_out:%m-%d %H:%M}, capped"
            else:
                label = f"OUT {t_out:%m-%d %H:%M}"
            ax.plot(t_min, y_fit, "-", lw=1.1, color=color, label=label)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(alpha=0.35)
        ax.legend(loc="best", fontsize=6, ncol=2)
    axs[0].set_title("Wet-bulb and humidity response after OUT events")
    axs[-1].set_xlabel("time since OUT (min)")
    fig.tight_layout()
    out_png = IMG_DIR / "Lab20260629_response_time_report.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"saved {out_png}")


def wet_bulb_settling(psy_df: pd.DataFrame) -> pd.DataFrame:
    df = psy_df.set_index("time")
    # Use only OUT-to-next-IN intervals for wet-bulb cooling; skip pre-fan windows.
    out_times = {ts for ts, kind, _label in events_as_timestamps() if kind == "OUT"}
    intervals = [
        (t0, t1, label, phase)
        for t0, t1, label, phase in air_windows()
        if t0 in out_times
    ]

    knee_frac = 0.10
    rows = []
    print(
        f"\n{'interval':<17s} {'N':>5s} {'Tw0':>6s} {'Tw_plat':>7s} "
        f"{'dTw':>6s} {'t_peak':>6s} {'t_knee':>6s} {'t63':>6s} "
        f"{'Td@plat':>8s} {'RH@plat':>8s}"
    )

    for ts0, ts1, label, phase in intervals:
        seg = df.loc[ts0:ts1]
        if len(seg) < 10:
            continue
        tw_seg = seg["T_wet_C"].to_numpy(dtype=float)
        td_seg = seg["T_dry_C"].to_numpy(dtype=float)
        rh_seg = seg["RH_pct"].to_numpy(dtype=float)
        t_rel = (seg.index - ts0).total_seconds().astype(float).to_numpy()
        tw_smooth = pd.Series(tw_seg).rolling(7, center=True, min_periods=1).mean().to_numpy()
        d_twdt = np.gradient(tw_smooth, t_rel)
        d_twdt = pd.Series(d_twdt).rolling(7, center=True, min_periods=1).mean().to_numpy()

        i_peak = int(np.argmin(d_twdt))
        peak_rate = d_twdt[i_peak]
        if peak_rate >= 0:
            t_peak = float(t_rel[i_peak])
            t_knee = None
            tw_plat = float(np.median(tw_smooth[-max(3, len(t_rel) // 10) :]))
        else:
            threshold = knee_frac * peak_rate
            after = np.where((t_rel > t_rel[i_peak]) & (d_twdt > threshold))[0]
            if after.size:
                i_knee = int(after[0])
                t_knee = float(t_rel[i_knee])
                tw_plat = float(np.median(tw_smooth[i_knee:]))
            else:
                t_knee = None
                tw_plat = float(np.median(tw_smooth[-max(3, len(t_rel) // 10) :]))
            t_peak = float(t_rel[i_peak])

        tw0 = float(tw_smooth[0])
        d_tw = tw0 - tw_plat
        target63 = tw0 - 0.63 * d_tw
        below63 = np.where(tw_smooth <= target63)[0]
        t63 = float(t_rel[below63[0]]) if below63.size and d_tw > 0 else None

        if t_knee is not None:
            plat_mask = t_rel >= t_knee
        else:
            plat_mask = np.zeros_like(t_rel, dtype=bool)
            plat_mask[-max(3, len(t_rel) // 10) :] = True
        td_at_plat = float(np.median(td_seg[plat_mask]))
        rh_at_plat = float(np.median(rh_seg[plat_mask]))
        t_settled_sec = t_knee if t_knee is not None else float(t_rel[-1])
        t_settled = ts0 + pd.Timedelta(seconds=t_settled_sec)

        def fmt(x: float | None) -> str:
            return "--" if x is None else f"{x:5.0f}s"

        print(
            f"{label:<17s} {len(seg):5d} {tw0:6.2f} {tw_plat:7.2f} "
            f"{d_tw:6.2f} {fmt(t_peak):>6s} {fmt(t_knee):>6s} {fmt(t63):>6s} "
            f"{td_at_plat:8.2f} {rh_at_plat:8.1f}"
        )

        rows.append(
            {
                "interval": label,
                "phase": phase,
                "t_out": ts0,
                "t_settled": t_settled,
                "Tw0": tw0,
                "Tw_plat": tw_plat,
                "dTw": d_tw,
                "t_peak_s": t_peak,
                "t_knee_s": t_knee,
                "t63_s": t63,
                "Td_at_plat": td_at_plat,
                "RH_at_plat": rh_at_plat,
            }
        )

    summary = pd.DataFrame(rows)
    out_csv = DATA_DIR / "wet_bulb_settling_20260629.csv"
    summary.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    wxt = load_wxt()
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    axes[0].plot(df.index, df["T_dry_C"], lw=0.8, color="C0", label="T_dry")
    axes[0].plot(df.index, df["T_wet_C"], lw=0.8, color="C1", label="T_wet")
    axes[0].plot(wxt.index, wxt["atmp"], lw=0.8, color="C3", label="WXT air T")
    axes[0].set_ylabel("T (deg C)")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].grid(alpha=0.4)
    axes[0].set_title("Lab20260629 wet-bulb cooling after OUT events")

    axes[1].plot(df.index, df["depression_C"], lw=0.8, color="C3")
    axes[1].set_ylabel("depression (deg C)")
    axes[1].grid(alpha=0.4)

    axes[2].plot(df.index, df["RH_pct"], lw=0.8, color="C0", label="psychrometric RH")
    axes[2].plot(wxt.index, wxt["hrh"], lw=0.8, color="C3", label="WXT520 RH")
    axes[2].set_ylabel("RH (%)")
    axes[2].set_ylim(0, 115)
    axes[2].legend(loc="best", fontsize=8)
    axes[2].grid(alpha=0.4)

    for row in rows:
        for ax in axes:
            ax.axvline(row["t_out"], color="tab:orange", ls="--", lw=0.6, alpha=0.7)
            ax.axvline(row["t_settled"], color="tab:green", ls=":", lw=0.6, alpha=0.7)
    for t0_s, _t1_s, _label, _phase in IN_WINDOWS:
        ts = pd.Timestamp(t0_s)
        for ax in axes:
            ax.axvline(ts, color="tab:blue", ls="--", lw=0.4, alpha=0.4)

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axes[-1].set_xlabel("time (UTC)")
    axes[-1].set_xlim(pd.Timestamp("2026-06-29 15:20"), pd.Timestamp("2026-06-30 16:45"))
    fig.tight_layout()
    out_png = IMG_DIR / "Lab20260629_temp_drops.png"
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f"saved {out_png}")
    return summary


def timeseries_overview(psy_df: pd.DataFrame, response_summary: pd.DataFrame) -> None:
    df = psy_df.set_index("time")
    wxt = load_wxt()
    wxt_pkl = pd.read_pickle(WXT_PKL)
    raw = wxt_pkl.loc[(wxt_pkl.index >= START) & (wxt_pkl.index <= END)]
    n_samples = len(raw["wspd11"].iloc[0])
    t_rep = np.repeat(raw.index.values, n_samples)
    wspd_raw = np.concatenate(raw["wspd11"].to_list())

    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
    axes[0].plot(df.index, df["T_dry_C"], lw=0.8, color="C0", label="RBR T_dry")
    axes[0].plot(df.index, df["T_wet_C"], lw=0.8, color="C1", label="RBR T_wet")
    axes[0].plot(wxt.index, wxt["atmp"], lw=0.8, color="C3", label="WXT air T")
    axes[0].set_ylabel("T (deg C)")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].grid(alpha=0.4)
    axes[0].set_title("RBRduo3 233860 + WXT520 ASWXT203 - Lab20260629")

    axes[1].plot(wxt.index, wxt["atmp"], color="C3", lw=0.8, label="atmp")
    axes[1].set_ylabel("WXT T (deg C)", color="C3")
    axes[1].tick_params(axis="y", labelcolor="C3")
    axes[1].grid(alpha=0.4)
    axb = axes[1].twinx()
    axb.plot(wxt.index, wxt["hrh"], color="C0", lw=0.8, label="hrh")
    axb.set_ylabel("WXT RH (%)", color="C0")
    axb.tick_params(axis="y", labelcolor="C0")
    h1, l1 = axes[1].get_legend_handles_labels()
    h2, l2 = axb.get_legend_handles_labels()
    axes[1].legend(h1 + h2, l1 + l2, loc="best", fontsize=8)

    axes[2].plot(wxt.index, wxt["bpr"], color="C2", lw=0.8, label="bpr")
    axes[2].set_ylabel("WXT pressure (mbar)")
    axes[2].legend(loc="best", fontsize=8)
    axes[2].grid(alpha=0.4)

    axes[3].plot(t_rep, wspd_raw, ".", ms=1.8, color="0.6", label="wspd11")
    axes[3].plot(wxt.index, wxt["wspd_min"], "-", lw=0.8, color="C2", label="wspd_min")
    axes[3].plot(wxt.index, wxt["wspd_max"], "-", lw=0.8, color="C3", label="wspd_max")
    axes[3].plot(wxt.index, wxt["wspd"], "-", lw=0.8, color="C0", label="wspd vector mean")
    axes[3].set_ylabel("wind speed (m/s)")
    axes[3].legend(loc="best", fontsize=8)
    axes[3].grid(alpha=0.4)

    mark_events(axes)
    axes[-1].set_xlim(START, END)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axes[-1].set_xlabel("time (UTC)")
    fig.tight_layout()
    out = DATA_DIR / "Lab20260629_RBR_WXT_timeseries.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"saved {out}")


def wxt_wind_plots() -> None:
    df = load_wxt()
    df_pkl = pd.read_pickle(WXT_PKL)
    df_pkl = df_pkl.loc[(df_pkl.index >= df.index.min()) & (df_pkl.index <= df.index.max())]

    n_samples = len(df_pkl["wdir11"].iloc[0])
    t_rep = np.repeat(df_pkl.index.values, n_samples)
    wdir = np.concatenate(df_pkl["wdir11"].to_list())
    wspd = np.concatenate(df_pkl["wspd11"].to_list())
    cmps = np.concatenate(df_pkl["compass11"].to_list())

    fig, axes = plt.subplots(5, 1, figsize=(11, 11), sharex=True)
    axes[0].plot(t_rep, wspd, ".", ms=2, color="0.6", label="wspd11")
    axes[0].plot(df.index, df["wspd"], "-", lw=0.8, color="C0", label="wspd vector mean")
    axes[0].set_ylabel("wind speed (m/s)")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].set_title(
        f"WXT520 / ASIMET PIC RevA203 - {WXT_DAT.name}\n"
        f"{df.index[0]} -> {df.index[-1]} UTC ({len(df)} records)"
    )

    axes[1].plot(t_rep, wdir, ".", ms=2, color="C0", label="wdir11")
    axes[1].set_ylabel("wind dir (deg)")
    axes[1].set_ylim(-10, 370)
    axes[1].legend(loc="best", fontsize=8)

    axes[2].plot(t_rep, cmps, ".", ms=2, color="C4", label="compass11")
    axes[2].plot(df.index, df["compass"], "-", lw=0.8, color="k", label="compass")
    axes[2].set_ylabel("compass (deg)")
    axes[2].set_ylim(-10, 370)
    axes[2].legend(loc="best", fontsize=8)

    axes[3].plot(df.index, df["atmp"], "-", color="C3", label="atmp (deg C)")
    axes[3].plot(df.index, df["hrh"], "-", color="C0", label="hrh (%)")
    axes[3].plot(df.index, df["bpr"] / 10.0, "-", color="C2", label="bpr/10")
    axes[3].set_ylabel("PTU")
    axes[3].legend(loc="best", fontsize=8)

    axes[4].plot(df.index, df["wndflag"], ".", color="C0", label="wndflag")
    axes[4].plot(df.index, df["rhtpflag"], ".", color="C3", label="rhtpflag")
    axes[4].plot(df.index, df["prcflag"], ".", color="C2", label="prcflag")
    axes[4].set_ylabel("status flags")
    axes[4].set_xlabel("time (UTC)")
    axes[4].legend(loc="best", fontsize=8)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    fig.tight_layout()
    out = IMG_DIR / "WXT_wind_samples_20260629.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"saved {out}")

    df2 = df.copy()
    df2["wspd_scalar_mean"] = df_pkl["wspd11"].apply(np.mean)
    df2["wspd11_min"] = df_pkl["wspd11"].apply(np.min)
    df2["wspd11_max"] = df_pkl["wspd11"].apply(np.max)

    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
    fig.suptitle(
        f"WXT520 winds - Lab 20260629\n"
        f"{df2.index[0]} -> {df2.index[-1]} UTC ({len(df2)} records)",
        fontsize=10,
    )

    axes[0].fill_between(df2.index, df2["wspd11_min"], df2["wspd11_max"], color="C0", alpha=0.25, label="wspd11 min/max")
    axes[0].plot(df2.index, df2["wspd_scalar_mean"], "-", lw=1.2, color="C0", label="wspd11 scalar mean")
    axes[0].plot(df2.index, df2["wspd"], "-", lw=1.0, color="C3", alpha=0.7, label="wspd vector mean")
    axes[0].set_ylabel("wind speed (m/s)")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].grid(alpha=0.4)

    axes[1].fill_between(df2.index, df2["wspd_min"], df2["wspd_max"], color="C2", alpha=0.25, label="instrument wspd_min/max")
    axes[1].plot(df2.index, df2["wspd_scalar_mean"], "-", lw=1.0, color="C0", alpha=0.7, label="wspd11 scalar mean")
    axes[1].set_ylabel("wind speed (m/s)")
    axes[1].legend(loc="best", fontsize=8)
    axes[1].grid(alpha=0.4)

    axes[2].scatter(df2.index, df2["wdir"], s=4, color="C1", label="wdir")
    axes[2].set_ylabel("wind dir (deg)")
    axes[2].set_ylim(-5, 365)
    axes[2].legend(loc="best", fontsize=8)
    axes[2].grid(alpha=0.4)

    axes[3].scatter(df2.index, df2["compass"], s=4, color="C4", label="compass")
    axes[3].set_ylabel("compass (deg)")
    axes[3].set_ylim(-5, 365)
    axes[3].set_xlabel("time (UTC)")
    axes[3].legend(loc="best", fontsize=8)
    axes[3].grid(alpha=0.4)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))

    fig.tight_layout()
    out = IMG_DIR / "WXT_wind_minmax_20260629.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"saved {out}")

    print("\n--- wind summary ---")
    for t0, t1, phase in ACTIVE_WINDOWS:
        sub = df2.loc[t0:t1]
        if sub.empty:
            continue
        print(
            f"{phase}: wspd11 scalar mean={sub['wspd_scalar_mean'].mean():.3f} m/s, "
            f"WXT RH={sub['hrh'].mean():.2f}%, WXT T={sub['atmp'].mean():.2f} deg C"
        )


def _latex_escape(text: str) -> str:
    return (
        str(text)
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("_", r"\_")
    )


def write_report(
    tank_summary: pd.DataFrame,
    wxt_comp: pd.DataFrame,
    response_summary: pd.DataFrame,
    settling_summary: pd.DataFrame,
) -> None:
    wxt = load_wxt()
    wxt_pkl = pd.read_pickle(WXT_PKL)
    wxt2 = wxt.copy()
    wxt2["wspd_scalar_mean"] = wxt_pkl["wspd11"].apply(np.mean)

    ambient_rows = []
    for t0, t1, phase in ACTIVE_WINDOWS:
        sub = wxt2.loc[t0:t1]
        ambient_rows.append(
            {
                "phase": phase,
                "T": sub["atmp"].mean(),
                "RH": sub["hrh"].mean(),
                "P": sub["bpr"].mean(),
                "wind": sub["wspd_scalar_mean"].mean(),
            }
        )

    valid = wxt_comp[wxt_comp["active_comparison_valid"]].copy()
    hum_rows = []
    for _t0, _t1, phase in ACTIVE_WINDOWS:
        sub = valid[valid["phase"] == phase]
        hum_rows.append(
            {
                "phase": phase,
                "N": len(sub),
                "wxt_rh": sub["hrh"].mean(),
                "psy_rh": sub["psy_RH_pct"].mean(),
                "d_rh": sub["dRH_psy_minus_wxt_pct"].mean(),
                "d_e": sub["de_psy_minus_wxt_hPa"].mean(),
                "d_q": sub["dspecific_humidity_psy_minus_wxt_gkg"].mean(),
            }
        )

    tau_rows = []
    for phase, sub in response_summary.groupby("phase"):
        tau_rows.append(
            {
                "phase": phase,
                "N": len(sub),
                "dry": sub["tau_T_dry_s"].median(),
                "wet": sub["tau_T_wet_s"].median(),
                "dep": sub["tau_depression_s"].median(),
                "rh": sub["tau_RH_psychro_s"].median(),
            }
        )

    settle_rows = []
    for phase, sub in settling_summary.groupby("phase"):
        settle_rows.append(
            {
                "phase": phase,
                "N": len(sub),
                "knee": sub["t_knee_s"].median(),
                "t63": sub["t63_s"].median(),
                "rh": sub["RH_at_plat"].median(),
            }
        )

    def ambient_table() -> str:
        lines = [
            r"\begin{tabular}{lrrrr}",
            r"\toprule",
            r"Phase & T (\si{\degreeCelsius}) & RH (\%) & p (mbar) & wind (m/s) \\",
            r"\midrule",
        ]
        for r in ambient_rows:
            lines.append(
                f"{_latex_escape(r['phase'])} & {r['T']:.2f} & {r['RH']:.1f} & "
                f"{r['P']:.2f} & {r['wind']:.2f} \\\\"
            )
        lines.extend([r"\bottomrule", r"\end{tabular}"])
        return "\n".join(lines)

    def tank_table() -> str:
        lines = [
            r"\begin{tabular}{llrrrr}",
            r"\toprule",
            r"Cycle & Phase & $T_\text{tank}$ & $T_d-T_\text{tank}$ & $T_w-T_\text{tank}$ & N \\",
            r"\midrule",
        ]
        for _, r in tank_summary.iterrows():
            lines.append(
                f"{_latex_escape(r['window'])} & {_latex_escape(str(r['phase']).replace('2026-06-', ''))} & "
                f"{r['T_tank_C']:.3f} & {r['dT_dry_minus_tank']:+.3f} & "
                f"{r['dT_wet_minus_tank']:+.3f} & {int(r['N'])} \\\\"
            )
        lines.extend([r"\bottomrule", r"\end{tabular}"])
        return "\n".join(lines)

    def humidity_table() -> str:
        lines = [
            r"\begin{tabular}{lrrrrr}",
            r"\toprule",
            r"Phase & N & WXT RH & psych RH & $\Delta$RH & $\Delta e$ \\",
            r"\midrule",
        ]
        for r in hum_rows:
            lines.append(
                f"{_latex_escape(r['phase'])} & {r['N']} & {r['wxt_rh']:.2f} & "
                f"{r['psy_rh']:.2f} & {r['d_rh']:+.2f} & {r['d_e']:+.2f} \\\\"
            )
        lines.extend([r"\bottomrule", r"\end{tabular}"])
        return "\n".join(lines)

    def tau_table() -> str:
        def tau_cell(value: float) -> str:
            if pd.isna(value):
                return "--"
            if value >= 0.98 * 7200.0:
                return r"$\geq7200^*$"
            return f"{value:.0f}"

        lines = [
            r"\begin{tabular}{lrrrrr}",
            r"\toprule",
            r"Phase & N & $\tau_{T_d}$ & $\tau_{T_w}$ & $\tau_\text{dep}$ & $\tau_\text{RH}$ \\",
            r"\midrule",
        ]
        for r in tau_rows:
            lines.append(
                f"{_latex_escape(r['phase'])} & {r['N']} & {tau_cell(r['dry'])} & "
                f"{tau_cell(r['wet'])} & {tau_cell(r['dep'])} & {tau_cell(r['rh'])} \\\\"
            )
        lines.extend([r"\bottomrule", r"\end{tabular}"])
        return "\n".join(lines)

    def settle_table() -> str:
        lines = [
            r"\begin{tabular}{lrrrr}",
            r"\toprule",
            r"Phase & N & median $t_\text{knee}$ (s) & median $t_{63}$ (s) & RH at plateau (\%) \\",
            r"\midrule",
        ]
        for r in settle_rows:
            knee = "--" if pd.isna(r["knee"]) else f"{r['knee']:.0f}"
            t63 = "--" if pd.isna(r["t63"]) else f"{r['t63']:.0f}"
            lines.append(
                f"{_latex_escape(r['phase'])} & {r['N']} & {knee} & {t63} & {r['rh']:.1f} \\\\"
            )
        lines.extend([r"\bottomrule", r"\end{tabular}"])
        return "\n".join(lines)

    tex = rf"""\documentclass[11pt,letterpaper]{{article}}

\usepackage[margin=1in]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{graphicx}}
\usepackage{{caption}}
\usepackage{{amsmath}}
\usepackage{{hyperref}}
\usepackage{{parskip}}
\usepackage{{microtype}}
\usepackage{{siunitx}}

\sisetup{{separate-uncertainty=true}}

\begin{{document}}

\noindent{{\LARGE\bfseries SHS Lab Test --- 2026-06-29\par}}
\vspace{{0.4ex}}
\noindent{{\large Clark South Ballast Tank, WHOI\par}}
\vspace{{1.0ex}}
\noindent Yu Gao; Emerson Hasbrouck; Tom Farrar; de Szoeke; Simon; Steven Jayne; Alexander Ekholm\par
\vspace{{0.4ex}}
\noindent Test date: 2026-06-29. The data record also includes the 2026-06-30 plastic-wick follow-up and post-test data through 2026-07-01 16:18 UTC.\par
\bigskip\hrule\bigskip

\section{{Setup}}

The test repeated the 2026-05-28 dunk-cycle workflow with two configurations.
On 2026-06-29 the SHS was run in the original wick setup; the fan was raised to
maximum after the first OUT event. The 2026-06-30 follow-up changed to a plastic
friction-fit wick at 12:55 UTC, with the fan turned on at the first IN.

\begin{{center}}
\begin{{tabular}}{{lll}}
\toprule
Instrument & ID & Rate \\
\midrule
RBRduo3 T.T (dry/wet bulb) & SN 233860 & 3 s \\
Vaisala WXT520 (met) & ASWXT203 & 1 min \\
Tank CT reference & --- & 1 Hz nominal \\
\bottomrule
\end{{tabular}}
\end{{center}}

\subsection*{{Event Log}}

\begin{{center}}
\begin{{tabular}}{{lllr}}
\toprule
Cycle & IN & OUT & Duration \\
\midrule
J29-IN-1 & 06-29 15:28 & 15:33 & 5 min \\
J29-IN-2 & 06-29 16:03 & 16:10 & 7 min \\
J29-IN-3 & 06-29 16:37 & 16:40 & 3 min \\
J29-IN-4 & 06-29 17:18 & 17:23 & 5 min \\
J29-IN-5 & 06-29 17:58 & 18:04 & 6 min \\
J29-IN-6 & 06-29 18:40 & 18:44 & 4 min \\
J29-IN-7 & 06-29 19:15 & 19:20 & 5 min \\
J29-IN-8 & 06-29 19:56 & 19:59 & 3 min \\
J30-IN-1 & 06-30 13:08 & 13:14 & 6 min \\
J30-IN-2 & 06-30 13:46 & 13:50 & 4 min \\
J30-IN-3 & 06-30 14:23 & 14:29 & 6 min \\
J30-IN-4 & 06-30 14:58 & 15:30 & 32 min \\
J30-IN-5 & 06-30 16:02 & 16:08 & 6 min \\
\bottomrule
\end{{tabular}}
\end{{center}}

\subsection*{{Ambient Conditions}}

\begin{{center}}
{ambient_table()}
\end{{center}}

\section{{Source Notes and Photos}}

\subsection*{{Testing Notes}}
The testing notes are reflected in the setup description and event log above.
They record RBR/WXT startup at 14:50 UTC, fan level 2 at 15:23 UTC, tank CT
logging at 15:27 UTC, the first SHS dunk at 15:28--15:33 UTC, and the fan
increase to maximum after the first OUT event with both probes facing the fan.
For the follow-up, the notes record the 12:55 UTC change to the plastic
friction-fit wick and the five subsequent dunk cycles.

\begin{{figure}}[p]
  \centering
  \includegraphics[width=0.48\linewidth,height=0.24\textheight,keepaspectratio]{{../../data/20260629Lab/IMG_7609.jpeg}}\hfill
  \includegraphics[width=0.48\linewidth,height=0.24\textheight,keepaspectratio]{{../../data/20260629Lab/IMG_7610.jpeg}}\\[0.8em]
  \includegraphics[width=0.31\linewidth,height=0.24\textheight,keepaspectratio]{{../../data/20260629Lab/IMG_7611.jpeg}}\hfill
  \includegraphics[width=0.31\linewidth,height=0.24\textheight,keepaspectratio]{{../../data/20260629Lab/IMG_7612.jpeg}}\hfill
  \includegraphics[width=0.31\linewidth,height=0.24\textheight,keepaspectratio]{{../../data/20260629Lab/IMG_7613.jpeg}}
  \caption{{Setup photographs from \texttt{{data/20260629Lab}}. Top: June 29
  original-wick layout with the WXT520, SHS wet/dry-bulb mast, fan, and ballast
  tank. Bottom: June 30 plastic-wick follow-up, including the friction-fit wick
  close-up and wider fan-facing instrument layout.}}
  \label{{fig:photo_log}}
\end{{figure}}

\section{{Humidity: WXT520 vs.\ Psychrometric}}

\begin{{figure}}[htbp]
  \centering
  \includegraphics[width=\linewidth]{{../../img/Lab20260629/Lab20260629_WXT_vs_psychrometric_humidity_report.png}}
  \caption{{Top: WXT520 RH and RBR wet/dry-bulb psychrometric RH. Bottom:
  psychrometric minus WXT RH difference and vapor-pressure difference. Grey
  bands mark in-water periods; the green dashed line marks the June 30
  plastic-wick change.}}
  \label{{fig:humidity}}
\end{{figure}}

\begin{{center}}
{humidity_table()}
\end{{center}}

The original-wick day tracks the WXT520 closely during active in-air windows.
The plastic-wick day shows a much larger psychrometric humidity offset, so that
configuration should be treated as unresolved until the wick geometry and wetting
state are checked. The full post-test record remains useful for seeing wick
drying and long-term drift, but the table above uses only active in-air samples
from the note-defined test windows.

\section{{Temperature and Tank Reference}}

\begin{{figure}}[htbp]
  \centering
  \includegraphics[width=\linewidth]{{../../img/Lab20260629/Lab20260629_SHS_vs_tank_report.png}}
  \caption{{Top: SHS dry- and wet-bulb temperatures with the in-tank CT
  reference. Bottom: SHS minus tank-temperature residuals during periods with CT
  coverage. Grey bands mark in-water periods; the green dashed line marks the
  June 30 plastic-wick change.}}
  \label{{fig:tank}}
\end{{figure}}

\begin{{center}}
{tank_table()}
\end{{center}}

\section{{Response Times}}

\begin{{figure}}[htbp]
  \centering
  \includegraphics[width=\linewidth]{{../../img/Lab20260629/Lab20260629_response_time_report.png}}
  \caption{{Wet-bulb temperature and psychrometric-RH trajectories aligned to
  each OUT event. Dots are observations and lines are exponential fits; labels
  list the OUT time and fitted time constant, or ``capped'' when the fit reached
  the upper bound. Fit windows run from OUT to
  $\min(\text{{OUT}}+25\,\text{{min}},\,\text{{next IN}})$.}}
  \label{{fig:response}}
\end{{figure}}

\begin{{center}}
{tau_table()}

\smallskip{{\footnotesize $^*$ Fit reached the 7200 s upper bound; treat those
values as non-exponential or unresolved rather than physical time constants.}}
\end{{center}}

\begin{{center}}
{settle_table()}
\end{{center}}

\section{{Diagnostics}}

The full multi-panel diagnostics are still generated for traceability, but are
omitted here to keep the report focused. They remain in
\texttt{{img/Lab20260629}} and \texttt{{data/20260629Lab}} alongside the
streamlined report figures.

\section{{Generated Files}}

\begin{{itemize}}
  \item \texttt{{data/20260629Lab/RBR\_SHS/psychrometrics\_20260629.csv}}
  \item \texttt{{data/20260629Lab/compare\_SHS\_tank\_20260629.csv}}
  \item \texttt{{data/20260629Lab/compare\_WXT\_psychrometric\_20260629.csv}}
  \item \texttt{{data/20260629Lab/response\_time\_20260629.csv}}
  \item \texttt{{data/20260629Lab/wet\_bulb\_settling\_20260629.csv}}
\end{{itemize}}

\end{{document}}
"""

    out_tex = DOC_DIR / "Lab20260629_report.tex"
    out_tex.write_text(tex)
    print(f"wrote {out_tex}")


def main() -> None:
    ensure_dirs()
    ensure_wxt_processed()
    ensure_rbr_nc()

    psy_df = psychrometrics()
    tank_summary = compare_shs_tank(psy_df)
    wxt_comp = compare_wxt_psychrometric(psy_df)
    response_summary = response_times(psy_df)
    settling_summary = wet_bulb_settling(psy_df)
    timeseries_overview(psy_df, response_summary)
    wxt_wind_plots()
    write_report(tank_summary, wxt_comp, response_summary, settling_summary)


if __name__ == "__main__":
    main()
