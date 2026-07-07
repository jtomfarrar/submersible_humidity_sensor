"""Python port of uop_decode.m (B. Greenwood / N. Galbraith, UOP/WHOI).

Orchestrates decoding of meteorological sensor data from the ASIMET instrument
suite. Run this script from the directory that contains the .DAT files, just
as you would run uop_decode.m from MATLAB.

Currently ported:
  - WXT520-SD  (Vaisala weather station, SD binary .DAT format)

Not yet ported (will warn and skip if listed):
  - Rotronic, Lascar, Campbell, SBE39, Airmar, ASIMet Logger, Standalone CF/SD

Usage
-----
Edit the USER CONFIG section below, cd to the data directory, then run::

    python /path/to/uop_decoder.py

Or pass a data directory explicitly::

    python /path/to/uop_decoder.py --datadir /path/to/WXT
"""
from __future__ import annotations

import argparse
import json
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# ── import WXT-SD decoder from the same package ───────────────────────────
import sys
sys.path.insert(0, str(Path(__file__).parent))
from decode_wxt520_sd import decode_wxt520

# ===========================================================================
#  USER CONFIG  (mirrors the "User defined" block in uop_decode.m)
# ===========================================================================
EXPERIMENT    = "SHS"
DEPLOYMENT    = "lab_2"
DUMPNAME      = "20260521_SHS-test_EGH"   # output sub-directory name
TSTART        = "2026-05-19"              # ISO start date/time (UTC, inclusive)
TEND          = "2026-05-21"              # ISO end   date/time (UTC, inclusive)
LATITUDE      =  41.52                    # deg N  (reserved for future use)
LONGITUDE     = -70.67                    # deg E  (reserved for future use)
AVERAGE_WINDS =  30                       # minutes to avg winds (0 = no avg)

# Sensor file stems (no extension).  Files must be in the data directory.
WXT_SD_FILES   = ["ASWXT102"]   # Vaisala WXT520-SD binary .DAT
ROTRONIC_FILES = []             # not yet ported
LASCAR_FILES   = []             # not yet ported
CAMPBELL_FILES = []             # not yet ported
SBE39_FILES    = []             # not yet ported
LOGGERS        = []             # ASIMet Logger .DAT – not yet ported
AIRMAR_DIR     = ""             # Airmar directory – not yet ported
# ===========================================================================


# ---------------------------------------------------------------------------
# Wind / direction averaging helpers (port of savg/davg from UOP toolbox)
# ---------------------------------------------------------------------------

def _scalar_avg(df: pd.DataFrame, col: str, rule: str) -> pd.Series:
    return df[col].resample(rule).mean()


def _vector_wind_avg(df: pd.DataFrame, rule: str) -> pd.DataFrame:
    """Vector-average wnde/wndn/wspd/wdir over *rule* (pandas offset string)."""
    ue = df["wnde"].resample(rule).mean()
    un = df["wndn"].resample(rule).mean()
    spd  = np.hypot(ue, un)
    dire = np.mod(np.degrees(np.arctan2(ue, un)), 360.0)
    return pd.DataFrame({"wnde": ue, "wndn": un, "wspd": spd, "wdir": dire})


def _directional_avg(series: pd.Series, rule: str) -> pd.Series:
    """Unit-vector mean for circular (directional) quantities."""
    rad = np.deg2rad(series)
    cu  = np.sin(rad).resample(rule).mean()
    cv  = np.cos(rad).resample(rule).mean()
    return np.mod(np.degrees(np.arctan2(cu, cv)), 360.0)


def _wind_avg(df: pd.DataFrame, minutes: int) -> pd.DataFrame:
    """Return a new DataFrame with winds and scalars averaged over *minutes*."""
    rule    = f"{minutes}min"
    dirlike = {"wdir", "wnde", "wndn", "wspd", "compass",
               "wdir11", "wspd11", "compass11"}

    scalar_cols = [c for c in df.columns if c not in dirlike]
    avg = df[scalar_cols].resample(rule).mean()

    vw = _vector_wind_avg(df, rule)
    for col in vw.columns:
        avg[col] = vw[col]

    avg["compass"] = _directional_avg(df["compass"], rule)
    return avg


# ---------------------------------------------------------------------------
# Per-sensor processor
# ---------------------------------------------------------------------------

def process_wxt_sd(
    infile: Path,
    outdir: Path,
    tstart: str,
    tend: str,
    average_winds: int,
) -> None:
    """Decode one WXT520-SD .DAT file and write outputs to *outdir*."""
    df, meta = decode_wxt520(infile)

    # trim to requested date range
    df = df.loc[tstart:tend]
    if df.empty:
        warnings.warn(f"{infile.name}: no records in [{tstart}, {tend}]")
        return

    meta["nrecs_decoded"] = int(len(df))
    meta["time_start"]    = df.index.min().isoformat()
    meta["time_end"]      = df.index.max().isoformat()

    scalar_cols = [c for c in df.columns
                   if c not in ("wdir11", "wspd11", "compass11")]

    # 1-minute CSV
    csv_path = outdir / f"{infile.stem}.csv"
    df[scalar_cols].to_csv(csv_path)

    # full pickle (preserves 11-sample arrays as Python lists)
    pkl_path = outdir / f"{infile.stem}.pkl"
    df.to_pickle(pkl_path)

    # optional wind-averaged CSV (mirrors MATLAB averaged_vars output)
    if average_winds:
        df_avg = _wind_avg(df, average_winds)
        avg_scalar = [c for c in df_avg.columns
                      if c not in ("wdir11", "wspd11", "compass11")]
        avg_csv = outdir / f"{infile.stem}_{average_winds}min_avg.csv"
        df_avg[avg_scalar].to_csv(avg_csv)
        print(f"  wrote {avg_csv}")

    # NetCDF (xarray optional)
    try:
        import xarray as xr
        ds = xr.Dataset.from_dataframe(df[scalar_cols])
        for col in ("wdir11", "wspd11", "compass11"):
            arr = np.vstack(df[col].to_numpy())
            ds[col] = (("time", "sample"), arr)
        for k, v in meta["instrument"].items():
            ds.attrs[f"instrument_{k}"] = v
        for k in ("firmware_version", "board_version", "module_sn",
                  "sensor_sn", "infile", "nrecs_decoded",
                  "bad_crc", "bad_time"):
            ds.attrs[k] = meta[k]
        nc_path = outdir / f"{infile.stem}.nc"
        ds.to_netcdf(nc_path)
        print(f"  wrote {nc_path}")
    except ImportError:
        pass

    meta_path = outdir / f"{infile.stem}_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2))

    print(f"  wrote {csv_path}")
    print(f"  wrote {pkl_path}")
    print(f"  wrote {meta_path}")
    print(f"  module_sn={meta['module_sn']}  sensor_sn={meta['sensor_sn']}  "
          f"fw={meta['firmware_version']}  "
          f"bad_crc={meta['bad_crc']}  bad_time={meta['bad_time']}")


# ---------------------------------------------------------------------------
# Main orchestration (mirrors uop_decode.m structure)
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--datadir", type=Path, default=None,
        help="directory containing .DAT files (default: current working dir)")
    ap.add_argument("--dumpname",  default=None,
                    help="output sub-directory name (overrides DUMPNAME)")
    ap.add_argument("--tstart",    default=None,
                    help="ISO start datetime UTC (overrides TSTART)")
    ap.add_argument("--tend",      default=None,
                    help="ISO end datetime UTC (overrides TEND)")
    ap.add_argument("--wxt",       default=None, nargs="+",
                    help="WXT file stem(s) without .DAT (overrides WXT_SD_FILES)")
    ap.add_argument("--avg-winds", type=int, default=None,
                    help="wind averaging window in minutes (overrides AVERAGE_WINDS)")
    args = ap.parse_args(argv)

    datadir       = args.datadir or Path.cwd()
    dumpname      = args.dumpname      or DUMPNAME
    tstart        = args.tstart        or TSTART
    tend          = args.tend          or TEND
    wxt_files     = args.wxt           or WXT_SD_FILES
    average_winds = args.avg_winds if args.avg_winds is not None else AVERAGE_WINDS

    # output layout mirrors uop_decode.m: DATADIR/processed/ and DATADIR/figs/
    outbase   = datadir / dumpname
    processed = outbase / "processed"
    figs      = outbase / "figs"
    for d in (outbase, processed, figs):
        d.mkdir(parents=True, exist_ok=True)

    # ── WXT520-SD ──────────────────────────────────────────────────────────
    for stem in wxt_files:
        infile = datadir / f"{stem}.DAT"
        if not infile.exists():
            warnings.warn(f"WXT file not found: {infile}")
            continue
        print(f"\n── WXT520-SD: {infile.name} ──")
        process_wxt_sd(infile, processed, tstart, tend, average_winds)
        try:
            shutil.copy2(infile, outbase / infile.name)
        except PermissionError:
            shutil.copy(infile, outbase / infile.name)  # skip metadata copy

    # ── Not-yet-ported sensors ─────────────────────────────────────────────
    for label, files in [
        ("Rotronic",      ROTRONIC_FILES),
        ("Lascar",        LASCAR_FILES),
        ("Campbell",      CAMPBELL_FILES),
        ("SBE39",         SBE39_FILES),
        ("ASIMet Logger", LOGGERS),
    ]:
        if files:
            warnings.warn(
                f"{label} decoder not yet ported to Python; skipping {files}")

    if AIRMAR_DIR:
        warnings.warn("Airmar decoder not yet ported to Python; skipping.")

    # copy this script to output directory (mirrors uop_decode.m self-archival)
    shutil.copy2(__file__, outbase / Path(__file__).name)

    print(f"\nDone. Output in: {outbase}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
