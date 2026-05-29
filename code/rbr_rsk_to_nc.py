"""Convert a simple two-channel RBR .rsk SQLite file to NetCDF."""
from __future__ import annotations

import argparse
from pathlib import Path

from rbr_rsk import write_rbr_netcdf


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rsk", type=Path, help="input .rsk file")
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="output NetCDF path (defaults to input path with .nc suffix)",
    )
    args = parser.parse_args(argv)

    out = args.out or args.rsk.with_suffix(".nc")
    written = write_rbr_netcdf(args.rsk, out)
    print(f"wrote {written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
