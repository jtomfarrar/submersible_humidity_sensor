"""Compute psychrometric humidity for the 2026-05-28 lab test from the RBR
wet/dry-bulb temperatures, and compare against the IN/OUT events recorded in
``data/20260528Lab/20260528_SHS-test-notes_EGH.txt``.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from rbr_rsk import write_rbr_netcdf


RD_OVER_RV = 287.04 / 461.5
A_FERREL = 0.000662  # Pa^-1 K^-1


def es(T_C: np.ndarray, p_Pa: float | np.ndarray) -> np.ndarray:
    """Saturation vapor pressure [Pa] over water (Buck 1981)."""
    p_hPa = p_Pa * 1e-2
    return 1e2 * 6.1121 * (1.0007 + 3.46e-8 * p_hPa) * np.exp(
        (17.502 * T_C) / (240.97 + T_C)
    )


def qs(p_Pa: float | np.ndarray, T_C: np.ndarray) -> np.ndarray:
    """Saturation specific humidity [kg/kg]."""
    e = es(T_C, p_Pa)
    return RD_OVER_RV * e / (p_Pa + (RD_OVER_RV - 1.0) * e)


def calc_psy(
    t_dry: np.ndarray, t_wet: np.ndarray, p_Pa: float | np.ndarray
) -> dict:
    """Return depression / RH / dew-point / specific-humidity arrays.

    ``p_Pa`` may be a scalar or an array matching ``t_dry``.
    """
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
    return dict(
        depression=depression,
        e=e,
        rh=rh,
        dew=dp,
        specific_humidity=specific_humidity,
    )


DAY = "2026-05-28"
NOTES = [
    ("13:08", "start/fan-on"),
    ("13:11", "IN"),
    ("13:27", "OUT"),
    ("13:56", "IN"),
    ("14:01", "OUT"),
    ("14:56", "IN"),
    ("15:01", "OUT"),
    ("15:42", "IN"),
    ("15:49", "OUT"),
    ("16:30", "IN"),
    ("16:34", "OUT"),
    ("18:30", "IN"),
    ("18:50", "OUT"),
    ("2026-05-29 11:10", "OFF"),
]


def note_time(value: str) -> pd.Timestamp:
    if value.startswith("20"):
        return pd.Timestamp(value)
    return pd.Timestamp(f"{DAY} {value}")


def state_after(kind: str, current: str) -> str:
    if "IN" in kind and "OUT" not in kind:
        return "water"
    if "OUT" in kind or "start" in kind:
        return "air"
    if kind == "OFF":
        return "off"
    return current


def main() -> None:
    here = Path(__file__).resolve().parent
    rbr_dir = here / "../data/20260528Lab/RBR"
    rbr_rsk = rbr_dir / "05_29_26_RBR.rsk"
    rbr_nc = rbr_dir / "05_29_26_RBR.nc"
    if not rbr_nc.exists():
        write_rbr_netcdf(rbr_rsk, rbr_nc)
        print(f"wrote {rbr_nc}")

    img_dir = here / "../img/Lab20260528"
    img_dir.mkdir(parents=True, exist_ok=True)

    ds = xr.open_dataset(rbr_nc)
    t_pd = pd.to_datetime(ds["time"].values)
    Td = ds["temperature"].values
    Tw = ds["temperature1"].values
    p_default = float(ds["default_atmospheric_pressure"].values[0]) * 1e4
    print(f"loaded {len(t_pd)} samples, default p_atm = {p_default / 100:.2f} hPa")

    wxt_csv = here / "../data/20260528Lab/WXT/processed/ASWXT203.csv"
    wxt = None
    if wxt_csv.exists():
        wxt = pd.read_csv(wxt_csv, parse_dates=["time"]).set_index("time")
        wxt = wxt.loc[(wxt.index >= t_pd[0]) & (wxt.index <= t_pd[-1])]

    # Use WXT-measured pressure (interpolated to RBR timestamps) when
    # available, else fall back to the RBR default scalar pressure.
    if wxt is not None and not wxt["bpr"].dropna().empty:
        wxt_p = wxt["bpr"].dropna()
        p_Pa = np.interp(
            mdates.date2num(t_pd),
            mdates.date2num(wxt_p.index.to_pydatetime()),
            wxt_p.values,
        ) * 100.0
        print(
            f"using WXT520 pressure for calc_psy: "
            f"min={p_Pa.min()/100:.2f} mean={p_Pa.mean()/100:.2f} "
            f"max={p_Pa.max()/100:.2f} hPa"
        )
    else:
        p_Pa = p_default
        print("using RBR default scalar pressure for calc_psy")

    psy = calc_psy(Td, Tw, p_Pa)
    events = [(note_time(h), kind) for h, kind in NOTES]

    in_air_mask = np.ones_like(Td, dtype=bool)
    in_air = True
    for ts, kind in events:
        if "IN" in kind and "OUT" not in kind:
            in_air = False
        elif "OUT" in kind or "start" in kind:
            in_air = True
        elif kind == "OFF":
            in_air = False
        in_air_mask[t_pd >= ts] = in_air

    print("\nMean psychrometric values for each interval defined by the notes:\n")
    print(
        f"{'interval':<21s} {'state':<6s} {'N':>5s} "
        f"{'Td':>6s} {'Tw':>6s} {'dep':>6s} "
        f"{'RH%':>5s} {'Tdew':>6s} {'spec_hum g/kg':>14s}"
    )

    t_last = pd.Timestamp(t_pd[-1])
    interval_start = pd.Timestamp(t_pd[0])
    state = "air"
    for ts, kind in events:
        if ts <= interval_start:
            state = state_after(kind, state)
            continue
        if ts > t_last:
            break
        sel = (t_pd >= interval_start) & (t_pd < ts)
        print_interval(interval_start, ts, state, sel, Td, Tw, psy)
        state = state_after(kind, state)
        interval_start = ts
    sel = t_pd >= interval_start
    print_interval(interval_start, t_last, state, sel, Td, Tw, psy)

    out_csv = rbr_dir / "psychrometrics_20260528.csv"
    df = pd.DataFrame(
        {
            "time": t_pd,
            "T_dry_C": Td,
            "T_wet_C": Tw,
            "depression_C": psy["depression"],
            "RH_pct": psy["rh"],
            "dew_point_C": psy["dew"],
            "specific_humidity_kg_per_kg": psy["specific_humidity"],
            "in_air": in_air_mask,
        }
    )
    df.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    fig, axs = plt.subplots(5, 1, figsize=(11, 13), sharex=True)

    axs[0].plot(t_pd, Td, label="T_dry (RBR temperature)", lw=0.8)
    axs[0].plot(t_pd, Tw, label="T_wet (RBR temperature1)", lw=0.8)
    if wxt is not None:
        axs[0].plot(wxt.index, wxt["atmp"], label="T_air (WXT520 atmp)",
                    lw=0.8, color="C3")
    axs[0].set_ylabel("T (deg C)")
    axs[0].legend(loc="best")
    axs[0].grid(alpha=0.4)
    axs[0].set_title(
        "SHS lab test 2026-05-28 - Ferrel psychrometric RH from RBRduo3 233860"
    )

    axs[1].plot(t_pd, psy["depression"], lw=0.8, color="C3")
    axs[1].set_ylabel("Td - Tw (deg C)")
    axs[1].grid(alpha=0.4)

    rh_air = np.where(in_air_mask, psy["rh"], np.nan)
    rh_sub = np.where(~in_air_mask, psy["rh"], np.nan)
    axs[2].plot(t_pd, rh_air, lw=0.8, color="C0", label="RH (in air, valid)")
    axs[2].plot(
        t_pd,
        rh_sub,
        lw=0.8,
        color="0.7",
        label="RH (sensor submerged - not physical)",
    )
    axs[2].set_ylabel("RH (%)")
    axs[2].set_ylim(0, 105)
    axs[2].grid(alpha=0.4)
    axs[2].legend(loc="best", fontsize=8)

    specific_humidity_air = np.where(
        in_air_mask, psy["specific_humidity"], np.nan
    ) * 1000.0
    axs[3].plot(t_pd, specific_humidity_air, lw=0.8, color="C2")
    axs[3].set_ylabel("specific humidity (g/kg, in air)")
    axs[3].grid(alpha=0.4)

    # Pressure used for the psychrometric calc -- WXT measured bpr (per-sample)
    # interpolated to RBR timestamps, with the RBR default scalar shown for
    # reference.
    p_used_hPa = (
        np.asarray(p_Pa) / 100.0 if np.ndim(p_Pa) else None
    )
    if p_used_hPa is not None:
        axs[4].plot(
            t_pd, p_used_hPa, lw=0.8, color="C0",
            label="pressure used in calc_psy (WXT bpr -> RBR times)",
        )
    axs[4].axhline(
        p_default / 100.0,
        color="C3",
        lw=1.0,
        ls="--",
        label=f"RBR default_atmospheric_pressure ({p_default / 100:.2f} hPa)",
    )
    if wxt is not None:
        axs[4].plot(
            wxt.index, wxt["bpr"], lw=0.8, color="C2", alpha=0.5,
            label="WXT520 measured pressure (bpr)",
        )
    axs[4].set_ylabel("pressure (hPa / mbar)")
    axs[4].legend(loc="best", fontsize=8)
    axs[4].grid(alpha=0.4)

    color_for = {
        "IN": "tab:blue",
        "OUT": "tab:orange",
        "start/fan-on": "k",
        "OFF": "k",
    }
    for ts, kind in events:
        c = color_for.get(kind, "gray")
        for ax in axs:
            ax.axvline(ts, color=c, ls="--", lw=0.6, alpha=0.6)
        if ts <= t_last:
            axs[0].annotate(
                kind,
                xy=(ts, axs[0].get_ylim()[1]),
                xytext=(0, -2),
                textcoords="offset points",
                rotation=90,
                ha="right",
                va="top",
                fontsize=7,
                color=c,
            )

    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axs[-1].set_xlabel("time (UTC)")
    fig.tight_layout()

    out_png = img_dir / "Lab20260528_SHS_humidity.png"
    fig.savefig(out_png, dpi=130)
    print(f"saved {out_png}")


def print_interval(
    ts0: pd.Timestamp,
    ts1: pd.Timestamp,
    state: str,
    sel: np.ndarray,
    Td: np.ndarray,
    Tw: np.ndarray,
    psy: dict,
) -> None:
    n = int(sel.sum())
    if n == 0:
        return
    label = f"{ts0.strftime('%m-%d %H:%M')}-{ts1.strftime('%m-%d %H:%M')}"
    print(
        f"{label:<21s} {state:<6s} {n:5d} "
        f"{np.nanmean(Td[sel]):6.2f} "
        f"{np.nanmean(Tw[sel]):6.2f} "
        f"{np.nanmean(psy['depression'][sel]):6.2f} "
        f"{np.nanmean(psy['rh'][sel]):5.1f} "
        f"{np.nanmean(psy['dew'][sel]):6.2f} "
        f"{np.nanmean(psy['specific_humidity'][sel]) * 1000:14.2f}"
    )


if __name__ == "__main__":
    main()
