import argparse

import numpy as np
import xarray as xr

from .. import codes
from .. import remo_ds as rds
from ..utils import encode
from . import _defaults as dflt
from .core import pressure_interpolation

plev_prec = np.float64

out_templates = {
    "input": "e{id}p_c{code:03d}_{date}.nc",
    "plev": "e{id}p_c{code:03d}_{plev:04d}_{date}.nc",
}

date_fmt = {"monthly": "%Y%m"}

INFO = """REMO pressure interpolaiton. Historically, REMO has its own interpolator although theoretically the algorithm is similar to what can be found, e.g., in metpy or GeoCAT. However, REMO’s algorithm includes a vertical spline fitting and requires temperature, pressure and orography as input."""


def generate_filename(id, time, code, plev=None):
    try:
        date = time.dt.strftime(date_fmt["monthly"]).values[0]
    except Exception:
        # brutal hack if time axis is absolute
        date = str(int(time.values[0]))[:6]
    if plev is not None:
        return out_templates["plev"].format(
            id=id, code=code, plev=int(plev.values / 100), date=date
        )
    else:
        return out_templates["input"].format(id=id, code=code, date=date)


def prsint_from_file(filename, varname, **args):
    ds = xr.open_dataset(filename)
    ds = rds.update_meta_info(ds)
    # ds = parse_dates(ds)
    result = pressure_interpolation(
        ds[varname],
        plev=args["plev"],
        t=ds.T,
        ps=ds.PS,
        orog=ds.FIB,
        a=ds.hyai,
        b=ds.hybi,
        keep_attrs=True,
    )
    return result


def write_output_variable(data_array, split, id):
    data_array = encode(data_array)
    if split is True:
        write_output_variable_one_per_plev(data_array, id)
    else:
        write_output_variable_like_input(data_array, id)


def write_output_variable_like_input(data_array, id):
    """writes data arrays to output files in the same structure as
    the input data.
    """
    time = data_array.time
    filename = generate_filename(
        id, time=time, code=codes.get_dict(data_array.name)["code"]
    )
    data_array.to_netcdf(filename)


def write_output_variable_one_per_plev(data_array, id):
    """writes data arrays to output files, one file per pressure level."""
    for plev in data_array.coords["plev"]:
        plev_data = data_array.sel(plev=[plev])
        time = data_array.time
        filename = generate_filename(
            id, time=time, code=codes.get_dict(data_array.name)["code"], plev=plev
        )
        # plev_data.to_netcdf(data_array.name+'_'+str(int(plev.values/100))+'.nc')
        plev_data.to_netcdf(filename)


def prsint(args):
    for f in args.input:
        for var in args.variables:
            output = prsint_from_file(
                f,
                varname=var,
                plev=args.plev,
            )
            write_output_variable(output, args.split, args.id)
    return 0


def prsint_parser():
    parser = argparse.ArgumentParser(description=INFO)
    parser.add_argument("input", metavar="input file", nargs="+")
    parser.add_argument(
        "-v",
        "--variables",
        dest="variables",
        nargs="+",
        help="list of variables to interpolate (default = {})".format(dflt.variables),
        default=dflt.variables,
    )
    parser.add_argument(
        "-p",
        "--plev",
        dest="plev",
        nargs="+",
        type=int,
        help="list of pressure levels to interpolate to (default = {})".format(
            dflt.plevs
        ),
        default=dflt.plevs,
    )
    parser.add_argument(
        "-s",
        "--split",
        dest="split",
        help="write one pressure level per output file",
        action="store_true",
    )
    parser.add_argument(
        "-id",
        "--id",
        dest="id",
        help="experiment id for output file naming",
        default="000000",
    )
    parser.add_argument(
        "-cdo",
        "--cdo_options",
        dest="cdo_options",
        help="options for using cdo to read input",
        default="",
    )
    return parser
