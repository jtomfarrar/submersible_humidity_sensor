"""Decode a Vaisala WXT520-SD raw binary .DAT file (UOP / ASIMET style 272-byte
records) into a pandas DataFrame / NetCDF / CSV.

This is a Python port of (the relevant parts of) get_sd_wxt520.m by
N. Galbraith / B. Greenwood (UOP, WHOI). It exists because the original MATLAB
script can fail in newer MATLABs at the `eval(savecmd)` line; the binary file
itself is fine.

Record layout (272 bytes, little-endian):
    bytes  0- 5 : 6 uchars  -> sec, min, hour, <unused>, day, month
                  (the MATLAB comment in get_sd_wxt520.m claims byte 3 is "day"
                  and byte 4 is "dow", but the actual datenum() call uses
                  time(5) as the day-of-month and ignores time(4), which is
                  what the raw bytes confirm.)
    bytes  6- 7 : uint16    -> year
    bytes  8-13 : 6 uchars  -> ASCII '272' (record size, ignored)
    bytes 14-15 : uint16    -> 272 (record byte size, ignored)
    bytes 16-207: 48 float32-> measurement payload (B array, 192 bytes)
    bytes 208-271: 64 uchars-> meta/firmware/flags (S array)
                              S[62:64] should be 0xA5 0xA5 (CRC marker)

B-array fields (48 float32):
     0-10  wdir11    (11 samples)
    11-21  wspd11    (11 samples)
    22     wspd_min
    23     wspd_max
    24-34  compass11 (11 samples)
    35     tilt_x
    36     tilt_y
    37     atmp
    38     hrh
    39     bpr
    40     precip
    41     rain_duration
    42     rain_intensity
    43     hail_accumulation
    44     hail_duration
    45     hail_intensity
    46     rain_peak_intensity
    47     hail_peak_intensity

S-array fields (64 uchar):
     0-19  version    (firmware string)
    20-35  brdversion (board firmware)
    36-39  modser     (module serial)
    40-47  sensor     (sensor serial)
    48-57  spare
    58     samp_count
    59     wndflag
    60     rhtpflag
    61     prcflag
    62-63  CRC marker (0xA5 0xA5)
"""

from __future__ import annotations

import argparse
import struct
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd


REC_SIZE = 272
CRC_BYTE = 0xA5


def _find_first_record(buf: bytes) -> int:
    """Scan for the 0xA5 0xA5 CRC marker, then back up REC_SIZE bytes to land
    on the start of the record whose end was just located."""
    for i in range(len(buf) - 1):
        if buf[i] == CRC_BYTE and buf[i + 1] == CRC_BYTE:
            # i+2 is just past the CRC. Step back REC_SIZE so cursor sits at
            # the start of the NEXT record (matches the MATLAB fseek(-272)).
            return (i + 2) - REC_SIZE
    raise RuntimeError("No 0xA5 0xA5 CRC marker found in file")


def _trim_cstr(b: bytes) -> str:
    """Strip trailing NULs / non-printable bytes from a fixed-width C string."""
    return b.split(b"\x00", 1)[0].decode("ascii", errors="replace").strip()


def decode_wxt520(infile: str | Path) -> tuple[pd.DataFrame, dict]:
    """Decode a WXT520-SD .DAT file.

    Returns
    -------
    df : pandas.DataFrame indexed by UTC time, with one row per 1-minute record.
    meta : dict with firmware / serial / file-level metadata.
    """
    infile = Path(infile)
    buf = infile.read_bytes()
    start = _find_first_record(buf)
    if start < 0:
        # very short file: first CRC is inside the first record, start at 0
        start = 0

    payload = buf[start:]
    nrecs = len(payload) // REC_SIZE
    if nrecs == 0:
        raise RuntimeError(f"File {infile} too small ({len(buf)} bytes)")

    print(f"Reading {infile.name}: {len(buf)} bytes, starting at offset "
          f"{start}, expecting up to {nrecs} records")

    records = []
    meta_fw = {"version": [], "brdversion": [], "modser": [], "sensor": []}
    bad_crc = 0
    bad_time = 0

    for i in range(nrecs):
        rec = payload[i * REC_SIZE:(i + 1) * REC_SIZE]
        if len(rec) < REC_SIZE:
            break

        # CRC check
        if rec[270] != CRC_BYTE or rec[271] != CRC_BYTE:
            bad_crc += 1
            continue

        sec, mn, hr, _unused, day, month = struct.unpack_from("<6B", rec, 0)
        (year,) = struct.unpack_from("<H", rec, 6)
        try:
            ts = datetime(year, month, day, hr, mn, sec)
        except ValueError:
            bad_time += 1
            continue

        B = np.frombuffer(rec, dtype="<f4", count=48, offset=16)
        S = np.frombuffer(rec, dtype=np.uint8, count=64, offset=208)

        rowd = {
            "time": ts,
            "wspd_min": float(B[22]),
            "wspd_max": float(B[23]),
            "tilt_x": float(B[35]),
            "tilt_y": float(B[36]),
            "atmp": float(B[37]),
            "hrh": float(B[38]),
            "bpr": float(B[39]),
            "precip": float(B[40]),
            "rain_duration": float(B[41]),
            "rain_intensity": float(B[42]),
            "hail_accumulation": float(B[43]),
            "hail_duration": float(B[44]),
            "hail_intensity": float(B[45]),
            "rain_peak_intensity": float(B[46]),
            "hail_peak_intensity": float(B[47]),
            "wndflag": int(S[59]),
            "rhtpflag": int(S[60]),
            "prcflag": int(S[61]),
        }
        wdir11 = B[0:11].astype(float)
        wspd11 = B[11:22].astype(float)
        compass11 = B[24:35].astype(float)

        # raw arrays as JSON-friendly lists (for traceability)
        rowd["wdir11"] = wdir11.tolist()
        rowd["wspd11"] = wspd11.tolist()
        rowd["compass11"] = compass11.tolist()

        # vector-average wind (meteorological convention, with compass rotation)
        tdir = np.mod(wdir11 + compass11, 360.0)
        spd_ok = wspd11 < 999.0
        if spd_ok.any():
            theta = np.deg2rad(tdir[spd_ok])
            spd = wspd11[spd_ok]
            wnde11 = spd * np.sin(theta)
            wndn11 = spd * np.cos(theta)
            wnde = float(np.mean(wnde11))
            wndn = float(np.mean(wndn11))
            wspd = float(np.hypot(wnde, wndn))
            wdir = float(np.mod(np.rad2deg(np.arctan2(wnde, wndn)), 360.0))
        else:
            wnde = wndn = wspd = wdir = np.nan

        # compass scalar average (unit-vector mean)
        cmp_rad = np.deg2rad(compass11)
        cu = float(np.mean(np.sin(cmp_rad)))
        cv = float(np.mean(np.cos(cmp_rad)))
        compass = float(np.mod(np.rad2deg(np.arctan2(cu, cv)), 360.0))

        rowd.update({"wnde": wnde, "wndn": wndn, "wspd": wspd,
                     "wdir": wdir, "compass": compass})
        records.append(rowd)

        meta_fw["version"].append(_trim_cstr(bytes(S[0:20])))
        meta_fw["brdversion"].append(_trim_cstr(bytes(S[20:36])))
        meta_fw["modser"].append(_trim_cstr(bytes(S[36:40])))
        meta_fw["sensor"].append(_trim_cstr(bytes(S[40:48])))

    if not records:
        raise RuntimeError("No valid records decoded")

    df = pd.DataFrame.from_records(records)
    df = df.sort_values("time").drop_duplicates(subset="time")
    # drop records whose timestamp goes backwards (matches MATLAB flag logic)
    keep = np.r_[True, np.diff(df["time"].values).astype("timedelta64[s]")
                 .astype(int) > 0]
    df = df.loc[keep].set_index("time")

    meta = {
        "infile": str(infile),
        "instrument": {
            "manufacturer": "Vaisala",
            "model": "WXT520",
            "model_vers": "SD",
        },
        "firmware_version": _mode(meta_fw["version"]),
        "board_version": _mode(meta_fw["brdversion"]),
        "module_sn": _mode(meta_fw["modser"]),
        "sensor_sn": _mode(meta_fw["sensor"]),
        "nrecs_decoded": int(len(df)),
        "bad_crc": int(bad_crc),
        "bad_time": int(bad_time),
        "time_start": df.index.min().isoformat(),
        "time_end": df.index.max().isoformat(),
    }
    return df, meta


def _mode(values: list[str]) -> str:
    if not values:
        return ""
    s = pd.Series(values)
    s = s[s != ""]
    if s.empty:
        return ""
    return str(s.mode().iloc[0])


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("infile", type=Path,
                   help="path to ASWXT*.DAT file")
    p.add_argument("--outdir", type=Path, default=None,
                   help="output directory (defaults to <infile_stem>_processed/)")
    p.add_argument("--tstart", type=str, default=None,
                   help="ISO start time to keep (UTC), inclusive")
    p.add_argument("--tend", type=str, default=None,
                   help="ISO end time to keep (UTC), inclusive")
    args = p.parse_args(argv)

    df, meta = decode_wxt520(args.infile)

    if args.tstart:
        df = df.loc[df.index >= pd.Timestamp(args.tstart)]
    if args.tend:
        df = df.loc[df.index <= pd.Timestamp(args.tend)]

    outdir = args.outdir or (args.infile.parent
                             / f"{args.infile.stem}_processed")
    outdir.mkdir(parents=True, exist_ok=True)

    # narrow CSV (drop the per-record 11-sample arrays which are awkward in CSV)
    scalar_cols = [c for c in df.columns
                   if c not in ("wdir11", "wspd11", "compass11")]
    csv_path = outdir / f"{args.infile.stem}.csv"
    df[scalar_cols].to_csv(csv_path)

    # full pickle preserves the 11-sample arrays as Python lists
    pkl_path = outdir / f"{args.infile.stem}.pkl"
    df.to_pickle(pkl_path)

    # try to also write a NetCDF (xarray is optional)
    try:
        import xarray as xr
        ds = xr.Dataset.from_dataframe(df[scalar_cols])
        # add the (time, 11) arrays as 2-D variables
        for col in ("wdir11", "wspd11", "compass11"):
            arr = np.vstack(df[col].to_numpy())
            ds[col] = (("time", "sample"), arr)
        for k, v in meta["instrument"].items():
            ds.attrs[f"instrument_{k}"] = v
        for k in ("firmware_version", "board_version", "module_sn",
                  "sensor_sn", "infile", "nrecs_decoded",
                  "bad_crc", "bad_time"):
            ds.attrs[k] = meta[k]
        nc_path = outdir / f"{args.infile.stem}.nc"
        ds.to_netcdf(nc_path)
        print(f"  wrote {nc_path}")
    except ImportError:
        print("  (xarray not available, skipping NetCDF output)")

    import json
    meta_path = outdir / f"{args.infile.stem}_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2))

    print(f"  wrote {csv_path}")
    print(f"  wrote {pkl_path}")
    print(f"  wrote {meta_path}")
    print(f"\nDecoded {len(df)} records "
          f"from {meta['time_start']} to {meta['time_end']}")
    print(f"  module_sn   = {meta['module_sn']}")
    print(f"  sensor_sn   = {meta['sensor_sn']}")
    print(f"  firmware    = {meta['firmware_version']}")
    print(f"  bad CRC     = {meta['bad_crc']}")
    print(f"  bad time    = {meta['bad_time']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
