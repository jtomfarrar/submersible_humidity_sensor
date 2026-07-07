"""Plot the 2026-05-28 RBR temperatures with WXT PTU, wind, and attitude."""
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import curve_fit

from decode_wxt520_sd import decode_wxt520
from rbr_rsk import write_rbr_netcdf
from vapo_sat_Lab20260528 import calc_psy as _calc_psy, NOTES, note_time


def _fit_tau(
    t_sec: np.ndarray, y: np.ndarray
) -> tuple[float, np.ndarray]:
    """Fit y = A + B·exp(−t/τ) to a step-response segment.
    Returns (tau_s, y_fitted).  t_sec[0] should be 0."""

    def _model(t, A, B, tau):
        return A + B * np.exp(-t / tau)

    A0 = float(y[-1])
    B0 = float(y[0] - y[-1])
    tau0 = max(float(t_sec[-1]) / 3.0, 10.0)
    try:
        popt, _ = curve_fit(
            _model, t_sec, y,
            p0=[A0, B0, tau0],
            bounds=([-np.inf, -np.inf, 1.0], [np.inf, np.inf, 7200.0]),
            maxfev=10000,
        )
        return float(popt[2]), _model(t_sec, *popt)
    except RuntimeError:
        return float("nan"), np.full(len(t_sec), np.nan)


REPO = Path(__file__).resolve().parents[1]
RBR_RSK = REPO / "data/20260528Lab/RBR/05_29_26_RBR.rsk"
RBR_NC = REPO / "data/20260528Lab/RBR/05_29_26_RBR.nc"
WXT_CSV = REPO / "data/20260528Lab/WXT/processed/ASWXT203.csv"
WXT_DAT = REPO / "data/20260528Lab/WXT/ASWXT203.DAT"
OUTPNG = REPO / "data/20260528Lab/Lab20260528_RBR_WXT_timeseries.png"
OUTPNG_TAU = REPO / "img/Lab20260528/Lab20260528_response_time.png"
TSTART = pd.Timestamp("2026-05-28 13:08")
TEND = pd.Timestamp("2026-05-29 11:10")


def main() -> None:
    if not RBR_NC.exists():
        write_rbr_netcdf(RBR_RSK, RBR_NC)
        print(f"wrote {RBR_NC}")

    ds = xr.open_dataset(RBR_NC)
    wxt = pd.read_csv(WXT_CSV, parse_dates=["time"]).set_index("time")
    wxt = wxt.loc[(wxt.index >= TSTART) & (wxt.index <= TEND)]

    # ── response-time estimation (all OUT events) ──────────────────────────
    # For each OUT event fit y = A + B·exp(−t/τ) over a 25-min window
    # (or until the next IN, whichever is shorter).
    FIT_WIN = pd.Timedelta(minutes=25)

    t_rbr = pd.DatetimeIndex(ds["time"].values)
    Td = ds["temperature"].values.astype(float)
    Tw = ds["temperature1"].values.astype(float)

    # Interpolate WXT pressure (mbar -> Pa) onto the RBR time grid
    p_Pa = np.interp(
        t_rbr.astype(np.int64),
        wxt.index.astype(np.int64),
        wxt["bpr"].values * 100.0,
    )
    psy = _calc_psy(Td, Tw, p_Pa)

    # Build (t_out, fit_end) pairs from NOTES
    out_events: list[tuple[pd.Timestamp, pd.Timestamp]] = []
    for i, (h, kind) in enumerate(NOTES):
        if kind == "OUT":
            t_out = note_time(h)
            t_next_in = None
            for hh, kk in NOTES[i + 1:]:
                if "IN" in kk and "OUT" not in kk:
                    t_next_in = note_time(hh)
                    break
            fit_end = t_out + FIT_WIN
            if t_next_in is not None and t_next_in < fit_end:
                fit_end = t_next_in
            out_events.append((t_out, fit_end))

    VAR_SPECS = [
        ("T_dry",      Td,                "T dry-bulb (deg C)"),
        ("T_wet",      Tw,                "T wet-bulb (deg C)"),
        ("depression", psy["depression"], "Depression T_d-T_w (deg C)"),
        ("RH_psychro", psy["rh"],         "Psychrometric RH (%)"),
    ]

    # all_tau_fits[lbl] = list of (t_arr, y_data, y_fit, tau) per OUT event
    all_tau_fits: dict[str, list] = {lbl: [] for lbl, *_ in VAR_SPECS}

    print("\n--- response-time estimates (all OUT events, 25-min window) ---")
    print(f"  {'event':<17s}  {'T_dry':>9s}  {'T_wet':>9s}  "
          f"{'depression':>11s}  {'RH_psychro':>11s}")
    for t_out, fit_end in out_events:
        fit_mask = (t_rbr >= t_out) & (t_rbr <= fit_end)
        t_sec = (t_rbr[fit_mask] - t_out).total_seconds().values.astype(float)
        row = []
        for lbl, arr, _ in VAR_SPECS:
            tau, y_fit = _fit_tau(t_sec, arr[fit_mask])
            all_tau_fits[lbl].append((t_rbr[fit_mask], arr[fit_mask], y_fit, tau))
            row.append(tau)
        print(
            f"  OUT {t_out.strftime('%H:%M')} "
            f"(win {int((fit_end - t_out).total_seconds() / 60)} min)  "
            + "  ".join(f"{tau:6.0f} s ({tau/60:4.1f} min)" for tau in row)
        )
    # ─────────────────────────────────────────────────────────────────────────

    # Convenience alias: first OUT event fits (used in timeseries overlay)
    tau_fits_first = {lbl: all_tau_fits[lbl][0] for lbl, *_ in VAR_SPECS}
    # ── response-time figure (all OUT events overlaid) ──────────────────────
    OUTPNG_TAU.parent.mkdir(parents=True, exist_ok=True)
    event_colors = [f"C{i}" for i in range(len(out_events))]
    fig_tau, axs_tau = plt.subplots(4, 1, figsize=(9, 11), sharex=True)
    for ax_t, (lbl, _, ylabel) in zip(axs_tau, VAR_SPECS):
        for ei, ((t_out, fit_end), color) in enumerate(
            zip(out_events, event_colors)
        ):
            t_arr, y_data, y_fit, tau = all_tau_fits[lbl][ei]
            t_min = (t_arr - t_out).total_seconds() / 60.0
            win_min = int((fit_end - t_out).total_seconds() / 60)
            ax_t.plot(t_min, y_data, ".", ms=3, color=color, alpha=0.6)
            lbl_str = (
                f"OUT {t_out.strftime('%H:%M')} ({win_min} min)  "
                f"τ={tau:.0f} s ({tau/60:.1f} min)"
            )
            ax_t.plot(t_min, y_fit, "-", lw=1.5, color=color, label=lbl_str)
        ax_t.set_ylabel(ylabel, fontsize=9)
        ax_t.legend(loc="best", fontsize=7)
        ax_t.grid(alpha=0.4)
    axs_tau[0].set_title(
        "RBRduo3 response after each OUT event — 2026-05-28 lab\n"
        "fit: $y = A + B\\,e^{-t/\\tau}$"
    )
    axs_tau[-1].set_xlabel("time since OUT (min)")
    fig_tau.tight_layout()
    fig_tau.savefig(OUTPNG_TAU, dpi=140)
    print(f"saved {OUTPNG_TAU}")
    # ─────────────────────────────────────────────────────────────────────────

    # Raw 11-sample/minute wind speeds straight from the binary record, so we
    # can show the per-minute spread on top of the WXT-reported min/max/avg.
    raw, _meta = decode_wxt520(WXT_DAT)
    raw = raw.loc[(raw.index >= TSTART) & (raw.index <= TEND)]
    n_samples = len(raw["wspd11"].iloc[0])
    t_rep = np.repeat(raw.index.values, n_samples)
    wspd_raw = np.concatenate(raw["wspd11"].to_list())

    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    ax = axes[0]
    ax.plot(ds["time"], ds["temperature"], lw=0.8, color="C0",
            label="RBR T_dry (temperature)")
    ax.plot(ds["time"], ds["temperature1"], lw=0.8, color="C1",
            label="RBR T_wet (temperature1)")
    ax.plot(wxt.index, wxt["atmp"], lw=0.8, color="C3",
            label="WXT520 air T (atmp)")
    # Overlay exponential fits for T_dry and T_wet (first OUT event)
    for lbl, color in [("T_dry", "C0"), ("T_wet", "C1")]:
        t_fit, _, y_fit, tau = tau_fits_first[lbl]
        ax.plot(t_fit, y_fit, "--", lw=1.5, color=color,
                label=f"{lbl} fit  τ={tau:.0f} s ({tau / 60:.1f} min)")
    ax.set_ylabel("T (deg C)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.4)
    ax.set_title("RBRduo3 233860 + WXT520 ASWXT203 - 2026-05-28 lab")

    ax = axes[1]
    ax.plot(wxt.index, wxt["atmp"], color="C3", lw=0.8, label="atmp")
    ax.set_ylabel("WXT air T (deg C)", color="C3")
    ax.tick_params(axis="y", labelcolor="C3")
    ax.grid(alpha=0.4)
    axb = ax.twinx()
    axb.plot(wxt.index, wxt["hrh"], color="C0", lw=0.8, label="hrh")
    axb.set_ylabel("WXT RH (%)", color="C0")
    axb.tick_params(axis="y", labelcolor="C0")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axb.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=8)

    ax = axes[2]
    ax.plot(wxt.index, wxt["bpr"], color="C2", lw=0.8, label="bpr")
    ax.set_ylabel("WXT pressure (mbar)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.4)

    ax = axes[3]
    # 11 raw 5-sec snapshots per minute, plus WXT-reported gust min/max/avg.
    # The snapshots are point samples while min/max/wspd are computed by the
    # WXT firmware from internal 3-sec running averages, so the snapshots
    # will not always fall between wspd_min and wspd_max.
    ax.plot(t_rep, wspd_raw, ".", ms=2, color="0.6",
            label="wspd11 (5-s snapshots, 11/min)")
    ax.plot(wxt.index, wxt["wspd_min"], "-", lw=0.8, color="C2",
            label="wspd_min (WXT 3-s gust min)")
    ax.plot(wxt.index, wxt["wspd_max"], "-", lw=0.8, color="C3",
            label="wspd_max (WXT 3-s gust max)")
    ax.plot(wxt.index, wxt["wspd"], "-", lw=0.8, color="C0",
            label="wspd (vector-averaged from wspd11)")
    ax.set_ylabel("wind speed (m/s)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.4)

    axes[-1].set_xlim(TSTART, TEND)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axes[-1].set_xlabel("time (UTC)")
    fig.tight_layout()
    fig.savefig(OUTPNG, dpi=130)
    print(f"saved {OUTPNG}")


if __name__ == "__main__":
    main()
