import argparse
from os import path as op
from pathlib import Path

import xarray as xr

soil_default = [
    "TSL",
    "TSN",
    "TD3",
    "TD4",
    "TD5",
    "GLAC",
    "TD",
    "TDCL",
    "WS",
    "WL",
    "SN",
]


def encode(ds):
    # we have to set the encoding here explicitly, otherwise xarray.to_netcdf will
    # encode missing values by NaN, which will crash REMO...
    fillvars = ["TSW", "SEAICE", "TSI"]
    for var in ds.data_vars:
        if var in fillvars:
            ds[var].encoding["_FillValue"] = 1.0e20
        else:
            ds[var].encoding["_FillValue"] = None
    for c in ds.coords:
        ds[c].encoding["_FillValue"] = None
    return ds


def get_output_filename(target):
    path = Path(target)
    return op.join(path.parent, path.stem + "_replaced" + path.suffix)


def replace_vars(target, source, vars, overwrite=False):
    tds = xr.open_dataset(target)
    sds = xr.open_dataset(source)
    tds = tds.merge(sds[vars], compat="override", join="override")
    if overwrite is True:
        fname = target
    else:
        fname = get_output_filename(target)
    print(f"writing to {fname}")
    encode(tds).to_netcdf(fname)


def replace_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("target", metavar="target", help="target file")
    parser.add_argument("source", metavar="source", help="source file")
    parser.add_argument(
        "-v",
        "--variables",
        dest="variables",
        nargs="+",
        help="list of variables to replace (default = {})".format(soil_default),
        default=soil_default,
    )
    parser.add_argument(
        "--overwrite",
        dest="overwrite",
        action="store_true",
        default=False,
        help="overwrite target file with new variables",
    )
    return parser


def replace_variables(args):
    replace_vars(args.target, args.source, args.variables, args.overwrite)
