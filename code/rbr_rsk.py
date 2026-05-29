"""Small RBR .rsk helpers for the lab-test analysis scripts.

The RBR ``.rsk`` files are SQLite databases.  These helpers read the measured
temperature channels directly, avoiding a dependency on ``pyrsktools`` for the
simple two-channel RBRduo3 lab-test files used here.
"""
from __future__ import annotations

import re
import sqlite3
from pathlib import Path

import pandas as pd
import xarray as xr


PARAMETER_NAMES = {
    "ATMOSPHERE": "default_atmospheric_pressure",
    "DENSITY": "default_density",
    "ALTITUDE": "default_altitude",
    "TEMPERATURE": "detault_temperature",
    "PRESSURE": "default_pressure",
    "CONDUCTIVITY": "default_conductivity",
    "SALINITY": "default_salinity",
    "SPECIFIC_CONDUCTIVITY_TEMPERATURE_CORRECTION": (
        "default_specific_conductivity_coefficient"
    ),
}


def _connect_readonly(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{path}?mode=ro", uri=True)


def _timestamp_ms_to_str(value: int | float | None) -> str:
    if value is None:
        return ""
    return pd.to_datetime(int(value), unit="ms").strftime("%Y-%m-%d %H:%M:%S")


def _safe_var_name(name: str) -> str:
    clean = re.sub(r"[^0-9A-Za-z_]+", "_", name.strip().lower()).strip("_")
    return clean or "channel"


def _channel_names(conn: sqlite3.Connection) -> list[str]:
    query = """
        select channels.longNamePlainText, instrumentChannels.channelOrder
        from instrumentChannels
        join channels on channels.channelID = instrumentChannels.channelID
        order by instrumentChannels.channelOrder
    """
    rows = conn.execute(query).fetchall()
    counts: dict[str, int] = {}
    names: list[str] = []
    for long_name, _order in rows:
        base = _safe_var_name(long_name)
        count = counts.get(base, 0)
        counts[base] = count + 1
        names.append(base if count == 0 else f"{base}{count}")
    return names


def read_rbr_rsk(path: str | Path) -> xr.Dataset:
    """Read a simple RBR ``.rsk`` SQLite file into an xarray Dataset."""
    path = Path(path)
    with _connect_readonly(path) as conn:
        data = pd.read_sql_query("select * from data order by tstamp", conn)
        if data.empty:
            raise RuntimeError(f"No samples found in {path}")

        times = pd.to_datetime(data.pop("tstamp"), unit="ms")
        names = _channel_names(conn)
        channel_cols = [c for c in data.columns if c.startswith("channel")]
        if len(names) != len(channel_cols):
            names = [f"channel{i + 1:02d}" for i in range(len(channel_cols))]

        variables = {
            name: ("time", data[col].to_numpy())
            for name, col in zip(names, channel_cols)
        }

        param_rows = conn.execute(
            """
            select parameterKeys.key, parameterKeys.value
            from parameterKeys
            join parameters
                on parameters.parameterID = parameterKeys.parameterID
            where parameters.parameterID = (
                select parameterID from parameters order by tstamp desc limit 1
            )
            """
        ).fetchall()
        params = {key: value for key, value in param_rows}
        for key, var_name in PARAMETER_NAMES.items():
            if key in params:
                try:
                    variables[var_name] = ("parameters", [float(params[key])])
                except ValueError:
                    pass

        instrument = conn.execute(
            """
            select serialID, model, firmwareVersion, firmwareType
            from instruments
            order by instrumentID
            limit 1
            """
        ).fetchone()
        deployment = conn.execute(
            """
            select timeOfDownload, name
            from deployments
            order by deploymentID
            limit 1
            """
        ).fetchone()
        ruskin = conn.execute(
            "select ruskinVersion from appSettings order by deploymentID limit 1"
        ).fetchone()
        db_info = conn.execute("select version from dbInfo limit 1").fetchone()

    ds = xr.Dataset(variables, coords={"time": times})
    serial, model, firmware_version, firmware_type = instrument or (
        "",
        "",
        "",
        "",
    )
    model_ascii = str(model).replace("³", "3")
    download_time = deployment[0] if deployment else None
    ds.attrs.update(
        {
            "Export Time": _timestamp_ms_to_str(download_time),
            "time coverage end": pd.Timestamp(times.max()).strftime(
                "%Y-%m-%d %H:%M:%S"
            ),
            "platform id": f"{model_ascii} {serial}".strip(),
            "Model": model_ascii,
            "Serial Number": str(serial),
            "Firmware Type": str(firmware_type),
            "Firmware Version": str(firmware_version),
            "Ruskin Version": ruskin[0] if ruskin else "",
            "File Version": db_info[0] if db_info else "",
            "standard name vocabulary": "CF Standard Name Table v65",
            "time coverage start": pd.Timestamp(times.min()).strftime(
                "%Y-%m-%d %H:%M:%S"
            ),
        }
    )
    return ds


def write_rbr_netcdf(rsk_path: str | Path, nc_path: str | Path) -> Path:
    """Decode ``rsk_path`` and write a NetCDF file at ``nc_path``."""
    nc_path = Path(nc_path)
    nc_path.parent.mkdir(parents=True, exist_ok=True)
    ds = read_rbr_rsk(rsk_path)
    ds.to_netcdf(nc_path)
    return nc_path
