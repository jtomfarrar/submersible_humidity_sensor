"""Plot WXT wind / compass / T-RH-P samples from the 2026-05-28 lab test."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from decode_wxt520_sd import decode_wxt520


REPO = Path(__file__).resolve().parents[1]
INFILE = REPO / "data/20260528Lab/WXT/ASWXT203.DAT"
OUTPNG = REPO / "img/Lab20260528/WXT_wind_samples_20260528.png"
TSTART = pd.Timestamp("2026-05-28 13:08")
TEND = pd.Timestamp("2026-05-29 11:10")
OUTPNG.parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    df, meta = decode_wxt520(INFILE)
    df = df.loc[(df.index >= TSTART) & (df.index <= TEND)]
    print(
        f"records={len(df)}  firmware={meta['firmware_version']}  "
        f"module_sn={meta['module_sn']}"
    )

    n_samples = len(df["wdir11"].iloc[0])
    t_rep = np.repeat(df.index.values, n_samples)
    wdir = np.concatenate(df["wdir11"].to_list())
    wspd = np.concatenate(df["wspd11"].to_list())
    cmps = np.concatenate(df["compass11"].to_list())

    fig, axes = plt.subplots(5, 1, figsize=(10, 11), sharex=True)

    ax = axes[0]
    ax.plot(t_rep, wspd, ".", ms=2, color="C0", label="wspd (11/min raw)")
    ax.plot(df.index, df["wspd_min"], "-", lw=0.8, color="C2", label="wspd_min")
    ax.plot(df.index, df["wspd_max"], "-", lw=0.8, color="C3", label="wspd_max")
    ax.set_ylabel("wind speed [m/s]")
    ax.legend(loc="best", fontsize=8)
    ax.set_title(
        f"WXT520 / ASIMET PIC RevA203 - {INFILE.name}\n"
        f"{df.index[0]} -> {df.index[-1]} UTC   "
        f"({len(df)} 1-min records, {n_samples} samples each)"
    )

    ax = axes[1]
    ax.plot(t_rep, wdir, ".", ms=2, color="C0", label="wdir (11/min raw)")
    ax.set_ylabel("wind dir [deg]")
    ax.set_ylim(-10, 370)
    ax.legend(loc="best", fontsize=8)

    ax = axes[2]
    ax.plot(t_rep, cmps, ".", ms=2, color="C4", label="compass (11/min raw)")
    ax.plot(df.index, df["compass"], "-", lw=0.8, color="k", label="compass (vec-avg)")
    ax.set_ylabel("compass [deg]")
    ax.set_ylim(-10, 370)
    ax.legend(loc="best", fontsize=8)

    ax = axes[3]
    ax.plot(df.index, df["atmp"], "-", color="C3", label="atmp [deg C]")
    ax.plot(df.index, df["hrh"], "-", color="C0", label="hrh [%]")
    ax.plot(df.index, df["bpr"] / 10.0, "-", color="C2", label="bpr/10 [mbar/10]")
    ax.set_ylabel("PTU")
    ax.legend(loc="best", fontsize=8)

    ax = axes[4]
    ax.plot(df.index, df["wndflag"], ".", color="C0", label="wndflag")
    ax.plot(df.index, df["rhtpflag"], ".", color="C3", label="rhtpflag")
    ax.plot(df.index, df["prcflag"], ".", color="C2", label="prcflag")
    ax.set_ylabel("status flags")
    ax.set_xlabel("time [UTC]")
    ax.legend(loc="best", fontsize=8)

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTPNG, dpi=140)
    print("wrote", OUTPNG)

    print("\n--- summary ---")
    print(
        "wspd  : min={:.3f}  mean={:.3f}  max={:.3f}  nonzero frac={:.3f}".format(
            wspd.min(), wspd.mean(), wspd.max(), (wspd != 0).mean()
        )
    )
    print(
        "wdir  : min={:.3f}  mean={:.3f}  max={:.3f}  nonzero frac={:.3f}".format(
            wdir.min(), wdir.mean(), wdir.max(), (wdir != 0).mean()
        )
    )
    print(
        "cmps  : min={:.3f}  mean={:.3f}  max={:.3f}  std={:.3f}".format(
            cmps.min(), cmps.mean(), cmps.max(), cmps.std()
        )
    )
    print(
        "atmp  : min={:.3f}  mean={:.3f}  max={:.3f}".format(
            df["atmp"].min(), df["atmp"].mean(), df["atmp"].max()
        )
    )
    print(
        "hrh   : min={:.3f}  mean={:.3f}  max={:.3f}".format(
            df["hrh"].min(), df["hrh"].mean(), df["hrh"].max()
        )
    )
    print(
        "bpr   : min={:.3f}  mean={:.3f}  max={:.3f}".format(
            df["bpr"].min(), df["bpr"].mean(), df["bpr"].max()
        )
    )
    print("wndflag value_counts:\n", df["wndflag"].value_counts().to_string())


if __name__ == "__main__":
    main()
