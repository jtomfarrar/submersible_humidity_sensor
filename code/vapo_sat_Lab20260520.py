"""Compute psychrometric humidity for the 2026-05-20 lab test from the RBR
wet/dry-bulb temperatures, and compare against the IN/OUT events recorded in
``data/20260520Lab/20260520_SHS-NOTES_EGH.txt``.

Uses the same Buck-1981 / Wexler saturation-vapor-pressure formulation as the
existing ``vapo_sat.py`` module and the same Ferrel-coefficient psychrometric
equation used in ``vapo_sat_Lab20250227.py``. The autograd import in the
original module is avoided by re-implementing ``es`` / ``qs`` here with plain
numpy (the autograd flavour is only needed for ``Twet_autodiff``).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


# ---------- constants (match vapo_sat.py) ---------------------------------
RD_OVER_RV = 287.04 / 461.5
A_FERREL = 0.000662  # Pa^-1 K^-1


def es(T_C: np.ndarray, p_Pa: float) -> np.ndarray:
    """Saturation vapor pressure [Pa] over water (Buck 1981, with enhancement)."""
    p_hPa = p_Pa * 1e-2
    return 1e2 * 6.1121 * (1.0007 + 3.46e-8 * p_hPa) * np.exp(
        (17.502 * T_C) / (240.97 + T_C)
    )


def qs(p_Pa: float, T_C: np.ndarray) -> np.ndarray:
    """Saturation specific humidity [kg/kg]."""
    e = es(T_C, p_Pa)
    return RD_OVER_RV * e / (p_Pa + (RD_OVER_RV - 1.0) * e)


def calc_psy(t_dry: np.ndarray, t_wet: np.ndarray, p_Pa: float) -> dict:
    """Return dict with depression / RH / dew-point / specific-humidity."""
    depression = t_dry - t_wet
    e_sat_dry = es(t_dry, p_Pa)
    e_sat_wet = es(t_wet, p_Pa)
    e = e_sat_wet - A_FERREL * p_Pa * depression          # actual vapor pressure
    rh = 100.0 * e / e_sat_dry
    q = qs(p_Pa, t_dry) * rh / 100.0
    # iterative dew-point (same scheme as vapo_sat_Lab20250227.py)
    dp = t_dry - 10.0
    for _ in range(40):
        e_dp = es(dp, p_Pa)
        e_dp1 = es(dp + 0.1, p_Pa)
        dp = dp + 0.2 * (e - e_dp) / (e_dp1 - e_dp) * 0.1
    return dict(depression=depression, e=e, rh=rh, dew=dp, q=q,
                e_sat_dry=e_sat_dry, e_sat_wet=e_sat_wet)


# ---------- IN/OUT events from 20260520_SHS-NOTES_EGH.txt ---------------
DAY = "2026-05-20"
NOTES = [
    ("14:15", "start"),         # RBR & WXT logging starts
    ("14:26", "IN"),            # SHS in water
    ("14:33", "fan-on/OUT"),    # fan on, sensor out of water
    ("15:07", "IN"),
    ("15:18", "OUT"),
    ("15:48", "IN"),
    ("15:53", "OUT"),
    ("17:02", "IN"),
    ("17:08", "OUT"),
    ("19:18", "IN"),
    ("19:25", "OUT"),
    ("20:59", "OFF"),
]


def main():
    here = Path(__file__).resolve().parent
    rbr_nc = here / "../data/20260520Lab/RBR/233860_20260521_1807.nc"
    img_dir = here / "../img"
    img_dir.mkdir(exist_ok=True)

    ds = xr.open_dataset(rbr_nc)
    t_pd = pd.to_datetime(ds["time"].values)
    Td = ds["temperature"].values        # dry bulb
    Tw = ds["temperature1"].values       # wet bulb
    p_atm = float(ds["default_atmospheric_pressure"].values[0]) * 1e4  # dbar -> Pa
    print(f"loaded {len(t_pd)} samples, p_atm = {p_atm/100:.2f} hPa")

    psy = calc_psy(Td, Tw, p_atm)

    # ---- build "in-air" mask from notes ---------------------------------
    events = [(pd.Timestamp(f"{DAY} {h}"), k) for h, k in NOTES]
    in_air_mask = np.ones_like(Td, dtype=bool)   # sensor in air before first IN
    state = True
    for ts, kind in events:
        idx = t_pd >= ts
        if "IN" in kind and "OUT" not in kind:
            state = False
        elif "OUT" in kind:
            state = True
        elif kind == "OFF":
            state = False
        in_air_mask[idx] = state

    # ---- per-interval summary -------------------------------------------
    print("\nMean psychrometric values for each interval defined by the notes:\n")
    print(f"{'interval':<21s} {'state':<6s} {'N':>5s} "
          f"{'Td':>6s} {'Tw':>6s} {'dep':>6s} "
          f"{'RH%':>5s} {'Tdew':>6s} {'q g/kg':>7s}")
    bounds = [t_pd[0]] + [ts for ts, _ in events] + [t_pd[-1]]
    state = "air"
    for i in range(len(bounds) - 1):
        ts0, ts1 = bounds[i], bounds[i + 1]
        sel = (t_pd >= ts0) & (t_pd < ts1)
        n = int(sel.sum())
        if n == 0:
            # advance state on the event that bounds this interval
            pass
        else:
            label = f"{ts0.strftime('%H:%M')}-{ts1.strftime('%H:%M')}"
            print(f"{label:<21s} {state:<6s} {n:5d} "
                  f"{np.nanmean(Td[sel]):6.2f} "
                  f"{np.nanmean(Tw[sel]):6.2f} "
                  f"{np.nanmean(psy['depression'][sel]):6.2f} "
                  f"{np.nanmean(psy['rh'][sel]):5.1f} "
                  f"{np.nanmean(psy['dew'][sel]):6.2f} "
                  f"{np.nanmean(psy['q'][sel])*1000:7.2f}")
        # update state for next interval based on event at ts1
        if i < len(events):
            kind = events[i][1] if i == 0 else events[i][1]
        # determine state for *next* interval from the event at ts1
        if i < len(events):
            ek = events[min(i, len(events) - 1)][1]
            if "IN" in ek and "OUT" not in ek:
                state = "water"
            elif "OUT" in ek:
                state = "air"
            elif ek == "OFF":
                state = "off"

    # ---- save CSV ---------------------------------------------------------
    out_csv = (rbr_nc.parent / "psychrometrics_20260520.csv")
    df = pd.DataFrame({
        "time": t_pd,
        "T_dry_C": Td,
        "T_wet_C": Tw,
        "depression_C": psy["depression"],
        "RH_pct": psy["rh"],
        "dew_point_C": psy["dew"],
        "q_kg_per_kg": psy["q"],
        "in_air": in_air_mask,
    })
    df.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    # ---- plot ------------------------------------------------------------
    fig, axs = plt.subplots(4, 1, figsize=(11, 11), sharex=True)

    axs[0].plot(t_pd, Td, label="T_dry (RBR temperature)", lw=0.8)
    axs[0].plot(t_pd, Tw, label="T_wet (RBR temperature1)", lw=0.8)
    axs[0].set_ylabel("T (°C)")
    axs[0].legend(loc="upper right")
    axs[0].grid(alpha=0.4)
    axs[0].set_title(
        "SHS lab test 2026-05-20 — Ferrel psychrometric RH from RBRduo3 233860"
    )

    axs[1].plot(t_pd, psy["depression"], lw=0.8, color="C3")
    axs[1].set_ylabel("Td − Tw  (°C)")
    axs[1].grid(alpha=0.4)

    rh_air = np.where(in_air_mask, psy["rh"], np.nan)
    rh_sub = np.where(~in_air_mask, psy["rh"], np.nan)
    axs[2].plot(t_pd, rh_air, lw=0.8, color="C0", label="RH (in air, valid)")
    axs[2].plot(t_pd, rh_sub, lw=0.8, color="0.7",
                label="RH (sensor submerged — not physical)")
    axs[2].set_ylabel("RH (%)")
    axs[2].set_ylim(0, 105)
    axs[2].grid(alpha=0.4)
    axs[2].legend(loc="lower right", fontsize=8)

    q_air = np.where(in_air_mask, psy["q"], np.nan) * 1000.0
    axs[3].plot(t_pd, q_air, lw=0.8, color="C2")
    axs[3].set_ylabel("q (g/kg, in air)")
    axs[3].grid(alpha=0.4)

    color_for = {"IN": "tab:blue", "OUT": "tab:orange",
                 "fan-on/OUT": "tab:red", "start": "k", "OFF": "k"}
    for ts, kind in events:
        c = color_for.get(kind, "gray")
        for ax in axs:
            ax.axvline(ts, color=c, ls="--", lw=0.6, alpha=0.6)
        axs[0].annotate(kind, xy=(ts, axs[0].get_ylim()[1]),
                        xytext=(0, -2), textcoords="offset points",
                        rotation=90, ha="right", va="top", fontsize=7, color=c)

    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[-1].set_xlabel("time (UTC)  2026-05-20")
    fig.tight_layout()

    out_png = img_dir / "Lab20260520_SHS_humidity.png"
    fig.savefig(out_png, dpi=130)
    print(f"saved {out_png}")


if __name__ == "__main__":
    main()
