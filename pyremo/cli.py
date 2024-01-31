"""Console script for pyremo."""

import argparse
import os
from os import path as op
from pathlib import Path

import xarray as xr

import pyremo as pr

from .core import cli as core_cli

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

fx_default = [
    "FIB",
    "BLA",
    "AZ0",
    "ALB",
    "VGRAT",
    "VAROR",
    "VLT",
    "FOREST",
    "FAO",
    "WSMX",
    "BETA",
    "WMINLOK",
    "WMAXLOK",
]


def create_parser():
    """Creates parser for command line tool."""
    return argparse.ArgumentParser()


def main():
    """Console script for pyremo."""
    parser = argparse.ArgumentParser()
    parser.add_argument("_", nargs="*")
    args = parser.parse_args()

    print("Arguments: " + str(args._))
    print("Replace this message by putting your code into pyremo.cli.main")
    return 0


def prsint():
    """Console script for prsint."""
    from .prsint import cli as prsint_cli

    parser = prsint_cli.prsint_parser()
    args = parser.parse_args()

    prsint_cli.prsint(args)
    return 0


def add_variables():
    """Console script for variable replacement."""
    parser = core_cli.replace_parser()
    args = parser.parse_args()

    core_cli.replace_variables(args)
    return 0


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


def replace_vars(target, source, vars, surflib=None, static=None, overwrite=False):
    tds = xr.open_dataset(target)
    sds = xr.open_dataset(source)
    tds = tds.merge(sds[vars], compat="override", join="override")
    if surflib is not None:
        surflib = pr.data.surflib(surflib)
        tds = tds.merge(surflib[static], compat="override", join="override")
    if overwrite is True:
        fname = target
    else:
        fname = get_output_filename(target)
    print(f"writing to {fname}")
    encode(tds).to_netcdf(fname)
    return fname


def replace_parser():
    parser = argparse.ArgumentParser(
        description="Replace variables in target files with data from source file. This can be used to, e.g., add soil variables from a REMO output file to a forcing file."
    )
    parser.add_argument("target", metavar="target", help="target file")
    parser.add_argument("source", metavar="source", help="source file")
    parser.add_argument(
        "-v",
        "--variables",
        dest="variables",
        nargs="+",
        help="list of variables to add or replace (default = {})".format(soil_default),
        default=soil_default,
    )
    parser.add_argument(
        "-s",
        "--surflib",
        dest="surflib",
        default=None,
        help="domain used for surface library",
    )
    parser.add_argument(
        "-fx",
        "--static",
        dest="static",
        nargs="+",
        help="list of variables to add or replace from soil library (default = {})".format(
            fx_default
        ),
        default=fx_default,
    )
    parser.add_argument(
        "-f",
        "--overwrite",
        dest="overwrite",
        action="store_true",
        default=False,
        help="overwrite target file with new variables",
    )
    return parser


def replace_variables(args):
    return replace_vars(
        args.target,
        args.source,
        args.variables,
        surflib=args.surflib,
        static=args.static,
        overwrite=args.overwrite,
    )


def main_analysis():
    import xarray as xr

    from .analysis import analysis, obs, plot
    from .analysis._parser import args

    time_range = slice(args.time_range[0], args.time_range[1])
    output_path = os.path.realpath(args.output_path)

    """
    Load Datasets. REMO and observational data
    """
    exp = analysis.RemoExperiment(args.input_path)
    remo_ds = exp.ds
    cru_ds = obs.cru_ts4(chunks={"time": 1, "lat": 360, "lon": 720})
    eobs_ds = obs.eobs(chunks={"time": 1, "latitude": 201, "longitude": 464})
    hyras_ds = obs.hyras()
    # hyras_ds.mask.plot(cmap='binary').savefig(os.path.join(output_path, 'HYRAS_mask.png'))

    """
    Convert data to into correct units and select time range
    """
    remo_tas = xr.merge(
        [remo_ds.TEMP2.sel(time=time_range) - 273.15, remo_ds.mask]
    ).rename({"TEMP2": "tas"})
    remo_tas.attrs["units"] = "Celsius"
    remo_pr = remo_ds.APRL.sel(time=time_range) + remo_ds.APRC.sel(time=time_range)
    remo_pr.name = "pr"
    remo_pr = xr.merge([remo_pr * 24, remo_ds.mask])

    eobs_tas = xr.merge([eobs_ds.tas.sel(time=time_range), eobs_ds.mask])
    eobs_pr = xr.merge([eobs_ds.pr.sel(time=time_range), eobs_ds.mask])

    cru_tas = xr.merge([cru_ds.tas.sel(time=time_range), cru_ds.mask])
    cru_pr = (
        cru_ds.pr.sel(time=time_range)
        / cru_ds.time.sel(time=time_range).dt.days_in_month
    )
    cru_pr.name = "pr"
    cru_pr.attrs["units"] = "mm/day"
    cru_pr = xr.merge([cru_pr, cru_ds.mask])

    hyras_tas = xr.merge([hyras_ds.tas.sel(time=time_range), hyras_ds.mask])
    hyras_pr = xr.merge([hyras_ds.pr.sel(time=time_range), hyras_ds.mask])

    """
    Plot seasonal differences
    """

    import cartopy.crs as ccrs

    pole = (
        remo_ds.rotated_latitude_longitude.grid_north_pole_longitude,
        remo_ds.rotated_latitude_longitude.grid_north_pole_latitude,
    )
    transform = ccrs.RotatedPole(*pole)

    # HYRAS
    extent = {
        "extents": [
            hyras_ds.lon.min(),
            hyras_ds.lon.max(),
            hyras_ds.lat.min(),
            hyras_ds.lat.max(),
        ]
    }

    # TAS
    compare_hyras_tas = analysis.compare_seasons(
        remo_tas, hyras_tas, regrid="ds2"
    ).compute()
    hyras_tas_season_plot = plot.plot_seasons(
        compare_hyras_tas.tas,
        transform=transform,
        extent=extent,
        borders=True,
        xlocs=range(-180, 180, 5),
        ylocs=range(-90, 90, 5),
        figsize=(14, 10),
        aspect="auto",
    )
    hyras_tas_season_plot.savefig(os.path.join(output_path, "TAS_REMO-HYRAS.png"))

    # PR
    compare_hyras_pr = analysis.compare_seasons(
        remo_pr, hyras_pr, regrid="ds2"
    ).compute()
    hyras_pr_season_plot = plot.plot_seasons(
        compare_hyras_pr.pr,
        transform=transform,
        extent=extent,
        borders=True,
        xlocs=range(-180, 180, 5),
        ylocs=range(-90, 90, 5),
        figsize=(14, 10),
        aspect="auto",
        vmax=5,
        vmin=-5,
        cmap="bwr_r",
    )
    hyras_pr_season_plot.savefig(os.path.join(output_path, "PR_REMO-HYRAS.png"))

    # CRU_TS4
    extent = {
        "extents": [
            remo_ds.rlon.min(),
            remo_ds.rlon.max(),
            remo_ds.rlat.min(),
            remo_ds.rlat.max(),
        ]
    }

    # TAS
    compare_cru_tas = analysis.compare_seasons(remo_tas, cru_tas).compute()
    cru_tas_season_plot = plot.plot_seasons(
        compare_cru_tas.tas,
        transform=transform,
        extent=extent,
        borders=True,
        xlocs=range(-180, 180, 10),
        ylocs=range(-90, 90, 10),
        figsize=(14, 10),
        aspect="auto",
    )
    cru_tas_season_plot.savefig(os.path.join(output_path, "TAS_REMO-CRU.png"))

    # PR
    compare_cru_pr = analysis.compare_seasons(remo_pr, cru_pr).compute()
    cru_pr_season_plot = plot.plot_seasons(
        compare_cru_pr.pr,
        transform=transform,
        extent=extent,
        borders=True,
        xlocs=range(-180, 180, 10),
        ylocs=range(-90, 90, 10),
        figsize=(14, 10),
        aspect="auto",
        vmin=-2,
        vmax=2,
        cmap="bwr_r",
    )

    cru_pr_season_plot.savefig(os.path.join(output_path, "PR_REMO-CRU.png"))

    # EOBS
    extent = {
        "extents": [
            remo_ds.rlon.min(),
            remo_ds.rlon.max(),
            remo_ds.rlat.min(),
            remo_ds.rlat.max(),
        ]
    }

    # TAS
    compare_eobs_tas = analysis.compare_seasons(remo_tas, eobs_tas).compute()
    eobs_tas_season_plot = plot.plot_seasons(
        compare_eobs_tas.tas,
        transform=transform,
        extent=extent,
        borders=True,
        xlocs=range(-180, 180, 10),
        ylocs=range(-90, 90, 10),
        figsize=(14, 10),
        aspect="auto",
    )
    eobs_tas_season_plot.savefig(os.path.join(output_path, "TAS_REMO-EOBS.png"))

    # PR
    compare_eobs_pr = analysis.compare_seasons(remo_pr, eobs_pr).compute()
    eobs_pr_season_plot = plot.plot_seasons(
        compare_eobs_pr.pr,
        transform=transform,
        extent=extent,
        borders=True,
        xlocs=range(-180, 180, 10),
        ylocs=range(-90, 90, 10),
        figsize=(14, 10),
        aspect="auto",
        vmin=-2,
        vmax=2,
        cmap="bwr_r",
    )

    eobs_pr_season_plot.savefig(os.path.join(output_path, "PR_REMO-EOBS.png"))
    return
