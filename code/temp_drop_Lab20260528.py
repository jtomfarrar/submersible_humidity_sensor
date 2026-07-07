"""Quantify how fast the SHS wet-bulb temperature equilibrates after each OUT
event during the 2026-05-28 lab test.

For every OUT->IN interval (the SHS is in air) we find the minimum wet-bulb
temperature reached and the time it took to get there, since the psychrometric
RH calculation is only meaningful once the wet bulb has finished cooling and is
on a quasi-steady plateau.

Produces a focused timeseries plot covering only the active in-air drops
(no overnight extension) and prints a table of the response times.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from rbr_rsk import write_rbr_netcdf
from vapo_sat_Lab20260528 import NOTES, calc_psy, note_time


REPO = Path(__file__).resolve().parents[1]
RBR_RSK = REPO / "data/20260528Lab/RBR/05_29_26_RBR.rsk"
RBR_NC = REPO / "data/20260528Lab/RBR/05_29_26_RBR.nc"
WXT_CSV = REPO / "data/20260528Lab/WXT/processed/ASWXT203.csv"
OUTPNG = REPO / "img/Lab20260528/Lab20260528_temp_drops.png"
OUTPNG.parent.mkdir(parents=True, exist_ok=True)

# Restrict to the active in-air period: from the first OUT to a bit after the
# last OUT. The overnight measurement (>20:00) is excluded.
TSTART = pd.Timestamp("2026-05-28 13:00")
TEND = pd.Timestamp("2026-05-28 19:30")


def out_in_intervals() -> list[tuple[pd.Timestamp, pd.Timestamp, str]]:
    """Return [(OUT time, next IN/OFF time, label), ...] within TSTART..TEND."""
    events = [(note_time(h), kind) for h, kind in NOTES]
    intervals = []
    for i, (ts, kind) in enumerate(events):
        if "OUT" not in kind and kind != "start/fan-on":
            continue
        # find next IN or OFF
        ts_end = None
        for ts2, k2 in events[i + 1:]:
            if "IN" in k2 and "OUT" not in k2:
                ts_end = ts2
                break
            if k2 == "OFF":
                ts_end = ts2
                break
        if ts_end is None:
            ts_end = TEND
        if ts > TEND or ts_end < TSTART:
            continue
        label = f"{ts.strftime('%H:%M')}-{ts_end.strftime('%H:%M')}"
        intervals.append((ts, min(ts_end, TEND), label))
    return intervals


def main() -> None:
    if not RBR_NC.exists():
        write_rbr_netcdf(RBR_RSK, RBR_NC)

    ds = xr.open_dataset(RBR_NC)
    t_pd = pd.to_datetime(ds["time"].values)
    Td = ds["temperature"].values
    Tw = ds["temperature1"].values
    p_atm = float(ds["default_atmospheric_pressure"].values[0]) * 1e4

    win = (t_pd >= TSTART) & (t_pd <= TEND)
    t_pd = t_pd[win]
    Td = Td[win]
    Tw = Tw[win]
    psy = calc_psy(Td, Tw, p_atm)

    wxt = None
    if WXT_CSV.exists():
        wxt = pd.read_csv(WXT_CSV, parse_dates=["time"]).set_index("time")
        wxt = wxt.loc[(wxt.index >= TSTART) & (wxt.index <= TEND)]

    intervals = out_in_intervals()

    # Per-OUT response statistics.
    # We want the inflection point / knee of the cooling curve: the time at
    # which the wet bulb transitions from active cooling into its quasi-steady
    # plateau.  This is more robust than picking the strict argmin (which is
    # noisy on flat plateaus) or a fixed fractional-amplitude crossing (which
    # is really a settling-time metric, not an inflection metric).
    #
    # Procedure per OUT->IN interval:
    #   1. Smooth Tw with a short rolling mean to suppress 1-sample noise.
    #   2. Compute dTw/dt (centered finite difference) and smooth it.
    #   3. Find the peak cooling rate (most-negative dTw/dt) and its time.
    #   4. The "knee" t_knee is the first time AFTER the peak cooling rate at
    #      which |dTw/dt| has decayed below 10 % of its peak magnitude --
    #      i.e. where the slope visibly flattens out.
    #   5. Plateau temperature Tw_plat = median of smoothed Tw from t_knee to
    #      the end of the interval.
    #   6. Also report t63 (e-folding) for reference.

    KNEE_FRAC = 0.10  # |dTw/dt| < 10% of peak  -> considered plateau

    print(
        f"{'interval':<13s} {'N':>5s} {'Tw0':>6s} {'Tw_plat':>7s} "
        f"{'dTw':>6s} {'t_peak':>6s} {'t_knee':>6s} {'t63':>6s} "
        f"{'Td@plat':>8s} {'RH@plat':>8s}"
    )
    summary_rows = []
    for ts0, ts1, label in intervals:
        sel = (t_pd >= ts0) & (t_pd <= ts1)
        if sel.sum() < 10:
            continue
        Tw_seg = Tw[sel]
        Td_seg = Td[sel]
        rh_seg = psy["rh"][sel]
        t_seg = t_pd[sel]
        t_rel = np.array(
            [(t - ts0).total_seconds() for t in t_seg], dtype=float
        )
        Tw_smooth = (
            pd.Series(Tw_seg).rolling(7, center=True, min_periods=1).mean()
            .to_numpy()
        )
        # centered finite difference, then smooth the derivative as well
        dTwdt_raw = np.gradient(Tw_smooth, t_rel)
        dTwdt = (
            pd.Series(dTwdt_raw).rolling(7, center=True, min_periods=1).mean()
            .to_numpy()
        )
        i_peak = int(np.argmin(dTwdt))  # most negative rate
        peak_rate = dTwdt[i_peak]
        if peak_rate >= 0:
            # never actually cooled -- skip the knee logic
            t_peak = t_rel[i_peak]
            t_knee = None
            Tw_plat = float(np.median(Tw_smooth[-max(3, len(t_rel) // 10):]))
        else:
            threshold = KNEE_FRAC * peak_rate  # less negative than peak
            after = np.where(
                (t_rel > t_rel[i_peak]) & (dTwdt > threshold)
            )[0]
            if after.size:
                i_knee = int(after[0])
                t_knee = float(t_rel[i_knee])
                Tw_plat = float(np.median(Tw_smooth[i_knee:]))
            else:
                t_knee = None
                Tw_plat = float(np.median(Tw_smooth[-max(3, len(t_rel) // 10):]))
            t_peak = float(t_rel[i_peak])

        Tw0 = float(Tw_smooth[0])
        dTw = Tw0 - Tw_plat

        target63 = Tw0 - 0.63 * dTw
        below63 = np.where(Tw_smooth <= target63)[0]
        t63 = float(t_rel[below63[0]]) if (below63.size and dTw > 0) else None

        # Td / RH on the plateau (use samples from t_knee onward, else last 10 %)
        if t_knee is not None:
            plat_mask = t_rel >= t_knee
        else:
            plat_mask = np.zeros_like(t_rel, dtype=bool)
            plat_mask[-max(3, len(t_rel) // 10):] = True
        Td_at_plat = float(np.median(Td_seg[plat_mask]))
        RH_at_plat = float(np.median(rh_seg[plat_mask]))

        def fmt(x: float | None) -> str:
            return "--" if x is None else f"{x:5.0f}s"

        print(
            f"{label:<13s} {sel.sum():5d} "
            f"{Tw0:6.2f} {Tw_plat:7.2f} {dTw:6.2f} "
            f"{fmt(t_peak):>6s} {fmt(t_knee):>6s} {fmt(t63):>6s} "
            f"{Td_at_plat:8.2f} {RH_at_plat:8.1f}"
        )

        t_settled_sec = t_knee if t_knee is not None else t_rel[-1]
        t_settled = ts0 + pd.Timedelta(seconds=t_settled_sec)
        summary_rows.append(
            dict(
                interval=label,
                t_out=ts0,
                t_settled=t_settled,
                Tw0=Tw0,
                Tw_plat=Tw_plat,
                dTw=dTw,
                t_peak=t_peak,
                t_knee=t_knee,
                t63=t63,
                Td_at_plat=Td_at_plat,
                RH_at_plat=RH_at_plat,
            )
        )

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)

    ax = axes[0]
    ax.plot(t_pd, Td, lw=0.8, color="C0", label="T_dry (RBR temperature)")
    ax.plot(t_pd, Tw, lw=0.8, color="C1", label="T_wet (RBR temperature1)")
    if wxt is not None:
        ax.plot(wxt.index, wxt["atmp"], lw=0.8, color="C3",
                label="T_air (WXT520 atmp)")
    ax.set_ylabel("T (deg C)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.4)
    ax.set_title(
        "SHS lab test 2026-05-28 - wet-bulb response after each OUT event "
        "(no overnight extension)"
    )

    ax = axes[1]
    ax.plot(t_pd, Td - Tw, lw=0.8, color="C3")
    ax.set_ylabel("depression Td - Tw (deg C)")
    ax.grid(alpha=0.4)

    ax = axes[2]
    ax.plot(t_pd, psy["rh"], lw=0.8, color="C0", label="psychrometric RH")
    if wxt is not None:
        ax.plot(wxt.index, wxt["hrh"], lw=0.8, color="C3",
                label="WXT520 RH (hrh)")
    ax.set_ylabel("RH (%)")
    ax.set_ylim(0, 110)
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.4)

    # Mark each OUT and the wet-bulb knee (inflection / end-of-cooling).
    for row in summary_rows:
        for a in axes:
            a.axvline(row["t_out"], color="tab:orange", ls="--", lw=0.6,
                      alpha=0.7)
            a.axvline(row["t_settled"], color="tab:green", ls=":", lw=0.6,
                      alpha=0.7)
        knee_label = (
            f"t_knee={row['t_knee']:.0f}s" if row["t_knee"] is not None
            else "no knee"
        )
        axes[0].annotate(
            knee_label,
            xy=(row["t_settled"], row["Tw_plat"]),
            xytext=(2, 6),
            textcoords="offset points",
            fontsize=7,
            color="tab:green",
        )
    # Mark IN events (sensor going back into water) as well, for context.
    for h, kind in NOTES:
        ts = note_time(h)
        if not (TSTART <= ts <= TEND):
            continue
        if "IN" in kind and "OUT" not in kind:
            for a in axes:
                a.axvline(ts, color="tab:blue", ls="--", lw=0.5, alpha=0.5)

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axes[-1].set_xlabel("time (UTC) on 2026-05-28")
    axes[-1].set_xlim(TSTART, TEND)
    fig.tight_layout()
    fig.savefig(OUTPNG, dpi=130)
    print(f"\nsaved {OUTPNG}")


if __name__ == "__main__":
    main()
