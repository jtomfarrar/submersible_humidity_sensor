"""Plot the 2026-05-28 RBR temperatures with WXT PTU and attitude channels."""
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

from rbr_rsk import write_rbr_netcdf


REPO = Path(__file__).resolve().parents[1]
RBR_RSK = REPO / "data/20260528Lab/RBR/05_29_26_RBR.rsk"
RBR_NC = REPO / "data/20260528Lab/RBR/05_29_26_RBR.nc"
WXT_CSV = REPO / "data/20260528Lab/WXT/processed/ASWXT203.csv"
OUTPNG = REPO / "data/20260528Lab/Lab20260528_RBR_WXT_timeseries.png"
TSTART = pd.Timestamp("2026-05-28 13:08")
TEND = pd.Timestamp("2026-05-29 11:10")


def main() -> None:
    if not RBR_NC.exists():
        write_rbr_netcdf(RBR_RSK, RBR_NC)
        print(f"wrote {RBR_NC}")

    ds = xr.open_dataset(RBR_NC)
    wxt = pd.read_csv(WXT_CSV, parse_dates=["time"]).set_index("time")

    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    ax = axes[0]
    ax.plot(ds["time"], ds["temperature"], lw=0.8, label="temperature")
    ax.plot(ds["time"], ds["temperature1"], lw=0.8, label="temperature1")
    ax.set_ylabel("RBR T (deg C)")
    ax.legend(loc="best")
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
    ax.plot(wxt.index, wxt["compass"], color="C4", lw=0.8, label="compass")
    ax.set_ylabel("WXT compass (deg)", color="C4")
    ax.tick_params(axis="y", labelcolor="C4")
    ax.grid(alpha=0.4)
    axb = ax.twinx()
    axb.plot(wxt.index, wxt["tilt_x"], color="C5", lw=0.8, label="tilt_x")
    axb.plot(wxt.index, wxt["tilt_y"], color="C6", lw=0.8, label="tilt_y")
    axb.set_ylabel("WXT tilt (deg)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axb.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=8)

    axes[-1].set_xlim(TSTART, TEND)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axes[-1].set_xlabel("time (UTC)")
    fig.tight_layout()
    fig.savefig(OUTPNG, dpi=130)
    print(f"saved {OUTPNG}")


if __name__ == "__main__":
    main()
