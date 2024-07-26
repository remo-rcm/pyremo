import glob
import os

# import subprocess
# import tempfile
from datetime import timedelta as td
from os import path as op
from warnings import warn

import cftime as cfdt
import dask
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from cdo import Cdo

from .core import (
    check_lev,
    convert_units,
    get_vc,
    horizontal_dims,
    map_sst,
    open_mfdataset,
)

cdo_exe = "cdo"

cdo_datetime_format = "%Y-%m-%dT%H:%M:%S"


@dask.delayed
def get_min_max_time(filename):
    with netCDF4.Dataset(filename) as ds:
        try:
            time_min = netCDF4.num2date(
                ds["time"][0], ds["time"].units, calendar=ds["time"].calendar
            )
            time_max = netCDF4.num2date(
                ds["time"][-1], ds["time"].units, calendar=ds["time"].calendar
            )
            return time_min, time_max
        except Exception:
            return np.nan, np.nan


def get_times_from_files(files):
    result = {}
    for f in files:
        result[f] = get_min_max_time(f)
    return result


def get_files(dirs):
    """get files from a list of directories"""
    if not isinstance(dirs, list):
        dirs = [dirs]
    files = []
    for d in dirs:
        d = op.abspath(d)
        # print('looking into', directory, op.isdir(directory), op.isfile(directory), Path(directory).is_dir())
        if op.isdir(d):
            files += glob.glob(os.path.join(d, "*.nc"))
        else:
            files += glob.glob(d)
    files.sort()
    return files


def create_var_df(var, data):
    df = pd.DataFrame.from_dict(data, orient="index", columns=["time_min", "time_max"])
    df.index.name = "path"
    df = df.reset_index()
    df["variable_id"] = var
    df.sort_values(by="time_min", inplace=True)
    return df


def create_df(data):
    return pd.concat([create_var_df(v, d) for v, d in data.items()], ignore_index=True)


def create_catalog(**args):
    """Scan directories for files and add time min and time max

    This function is supposed to scan directories for global model input files.
    It will create a catalog that contains the path to the file and the
    min and max times of that file.

    """
    files = {}
    print(args)
    for v, d in args.items():
        files[v] = get_files(d)
    for v, fs in files.items():
        files[v] = get_times_from_files(fs)
    return files


def to_datetime(time):
    try:
        return pd.to_datetime(str(int(time)))
    except Exception:
        return time


def to_cfdatetime(time, calendar="standard"):
    try:
        return xr.cftime_range(start=time, end=time, calendar=calendar)[0]
    except Exception:
        return np.nan


def search_df(df, **kwargs):
    """Search dataframe by arbitray conditions

    Converts kwargs to pandas search conditions. If kwargs is a list,
    pandas isin is used as condition.

    """
    condition_list = []
    for key, item in kwargs.items():
        if isinstance(item, list):
            cond = "(df['{0}'].isin({1}))".format(key, repr(item))
        else:
            cond = "(df['{0}'] == {1})".format(key, repr(item))
        condition_list.append(cond)
    conditions = " & ".join(condition_list)
    return df[eval(conditions)]


def get_var_by_time(df, datetime=None, **kwargs):
    df = search_df(df, **kwargs)
    if datetime:
        df = df[(datetime >= df.time_min) & (datetime <= df.time_max)]
    return df


def get_sst_times(date):
    """Get daily dates from which the SST is interpolated in time"""
    if date.hour == 12:
        # no interpolation neccessary
        return (date,)
    # if dt.hour > 12:
    #     # interpolate between today and tomorrow
    #     dt1 = dt + td(days=1)
    # else:
    #     # interpolated betweend today and yesterday
    #     dt1 = dt
    #     dt = dt + td(days=-1)
    # return (
    #     cfdt.datetime(dt.year, dt.month, dt.day, 12, calendar=cal),
    #     cfdt.datetime(dt1.year, dt1.month, dt1.day, 12, calendar=cal),
    # )

    days = (date + td(days=i) for i in range(-2, 3))
    return [
        cfdt.datetime(d.year, d.month, d.day, 12, calendar="standard") for d in days
    ]


def cdo_call(self, options="", op="", input="", output="temp", print_command=True):
    pass


class CFModelSelector:
    """CF Model Selector

    The CF model selector class simply selects files by a certain timestep.

    """

    def __init__(self, df, calendar="standard", **kwargs):
        df = df.copy()
        if kwargs:
            df = search_df(df, **kwargs)
        self.calendar = calendar
        self.tempfiles = []
        self.df = self._update_time(df)

    def _update_time(self, df):
        df["time_min"] = df["time_min"].apply(to_cfdatetime, calendar=self.calendar)
        df["time_max"] = df["time_max"].apply(to_cfdatetime, calendar=self.calendar)
        return df

    def get_file(self, datetime=None, **kwargs):
        sel = get_var_by_time(self.df, datetime=datetime, **kwargs)
        if len(sel.index) > 1:
            return list(sel.path)
        if sel.empty:
            raise FileNotFoundError(
                "no file found: {}, datetime: {}".format(kwargs, datetime)
            )
        return sel.iloc[0].path


def gfile(ds, ref_ds=None, tos=None, attrs=None, use_cftime=True, invertlev=None):
    """Creates a global dataset ready for preprocessing.

    This function creates a homogenized global dataset. If neccessary,
    units are converted and the sea surface temperature ``tos`` is
    interpolated spatially and temporally to the atmospheric grid.

    Parameters
    ----------
    ds : xarray.Dataset or dict of filenames
        Input dataset from a global model according to CF conventions.

    ref_ds : xarray.Dataset
        Reference datasets that is used for determining the grid and vertical
        coordinates and the global attributes. If ``ref_ds=None``, ``ta`` from
        the input dataset is used as a reference.

    attrs:
        Global attributes for the output dataset. If ``attrs=None``, the global
        attributes from ``ref_ds`` are used.

    Returns
    -------
    gfile : xarray.Dataset
        Global dataset ready for preprocessing.

    """
    if isinstance(ds, dict):
        ds = open_datasets(ds, ref_ds)
    else:
        ds = ds.copy()
        ds["akgm"], ds["bkgm"] = get_vc(ds, invertlev)
        ds = check_lev(ds, invertlev)

    if tos is not None:
        ds["tos"] = map_sst(tos, ds)
    # ensure correct units
    ds = convert_units(ds)
    if "sftlf" in ds:
        ds["sftlf"] = np.around(ds.sftlf)
    if attrs is None:
        attrs = ds.attrs
    if use_cftime is True:
        # what for https://github.com/pydata/xarray/pull/7399
        # return ds.convert_calendar(ds.time.dt.calendar, use_cftime=True)
        return ds.assign_coords(
            time=ds.time.convert_calendar(ds.time.dt.calendar, use_cftime=True).time
        )
    return ds


class GFile:
    dynamics = ["ta", "ua", "va", "ps", "hus"]
    fx = ["orog", "sftlf"]
    sst = "tos"
    all_vars = dynamics + fx  # + sst

    def __init__(self, df=None, calendar=None, scratch=None, **kwargs):
        if calendar is None:
            try:
                self.calendar = df.time_min.iloc[0].calendar
            except Exception:
                self.calendar = "standard"
        self.selector = CFModelSelector(df=df, calendar=self.calendar, **kwargs)
        self.scratch = scratch
        if self.scratch is None:
            try:
                self.scratch = os.path.join(os.environ["SCRATCH"], ".cf-selector")
            except Exception:
                pass

        self.tos_regridder = None

    def get_files(self, variables, datetime=None, **kwargs):
        if datetime:
            datetime = to_cfdatetime(datetime, self.calendar)
        files = {
            var: self.selector.get_file(variable_id=var, datetime=datetime, **kwargs)
            for var in variables
        }
        return files

    def extract_timestep(self, datetime=None, **kwargs):
        pass

    def extract_dynamic_timesteps(self, datetime=None, **kwargs):
        """Extract timesteps for a certain date from input files."""

        datetime = to_cfdatetime(datetime, self.calendar)
        files = self.get_files(self.dynamics, datetime=datetime, **kwargs)
        if datetime is None:
            return files
        if self.scratch is not None:
            cdo = Cdo(tempdir=self.scratch)
        else:
            cdo = Cdo()
        return {
            v: cdo.seldate(
                datetime.strftime(cdo_datetime_format), input=f, returnXDataset=True
            ).load()
            for v, f in files.items()
        }

    def extract_data(self, datetime, **kwargs):
        files = self.extract_dynamic_timesteps(datetime, **kwargs)
        files.update(
            {
                var: xr.open_dataset(f).load()
                for var, f in self.get_files(self.fx, datetime=None, **kwargs).items()
            }
        )

        return xr.merge(files.values(), compat="override", join="override")

    def get_sst(self, datetime, atmo_grid=None):
        """Extract and interpolate SST in space and time."""

        datetime = to_cfdatetime(datetime, self.calendar)
        times = get_sst_times(datetime)
        # files = {
        #    t: self.selector.get_file(variable_id=self.sst, datetime=t) for t in times
        # }
        files = {}
        for t in times:
            print(t)
            try:
                f = self.selector.get_file(variable_id=self.sst, datetime=t)
                files[t] = f
                print(f)
            except FileNotFoundError:
                warn(f"sst not found for {datetime}, will extrapolate...")
                pass

        if self.scratch is not None:
            cdo = Cdo(tempdir=self.scratch)
        else:
            cdo = Cdo()
        sst_extract = [
            cdo.seldate(
                t.strftime(cdo_datetime_format), input=f, returnXDataset=True
            ).load()
            for t, f in files.items()
        ]
        sst_da = xr.merge(sst_extract)[
            self.sst
        ]  # xr.open_mfdataset(sst_extract, use_cftime=True)
        if len(sst_da.time) > 1:
            sst_da = sst_da.interp(
                time=datetime.strftime(cdo_datetime_format),
                method="linear",
                kwargs={"fill_value": "extrapolate"},
            )
        if atmo_grid is None:
            return sst_da
        return self.regrid_to_atmosphere(sst_da.squeeze(drop=True), atmo_grid)

    def regrid_to_atmosphere(self, da, atmo_grid):
        """Regrid tos to atmoshperic grid."""
        import xesmf as xe

        attrs = da.attrs
        atmo_grid = atmo_grid.copy()
        atmo_grid["mask"] = ~(atmo_grid.sftlf > 0).squeeze(drop=True)

        ds = da.to_dataset()

        if not self.tos_regridder:
            ds["mask"] = ~ds.tos.isnull().squeeze(drop=True)
            print("creating tos regridder")
            self.tos_regridder = xe.Regridder(
                ds, atmo_grid, method="nearest_s2d", extrap_method="nearest_s2d"
            )

        out = self.tos_regridder(da)
        out.attrs = attrs
        return out

    def gfile(self, datetime, sst=True, **kwargs):
        """Creates a gfile from CF input data."""
        print(f"extracting: {datetime}")
        gds = self.extract_data(datetime=datetime, **kwargs)
        if sst is True:
            print("interpolating SST")
            gds[self.sst] = self.get_sst(datetime, gds)
        if "variable_id" in gds.attrs:
            del gds.attrs["variable_id"]
        return gfile(gds)


def open_datasets(datasets, ref_ds=None):
    """Creates a virtual gfile"""
    if ref_ds is None:
        try:
            ref_ds = open_mfdataset(datasets["ta"])
        except Exception:
            raise Exception("ta is required in the datasets dict if no ref_ds is given")
    lon, lat = horizontal_dims(ref_ds)
    # ak_bnds, bk_bnds = get_ab_bnds(ref_ds)
    dsets = []
    for var, f in datasets.items():
        try:
            da = open_mfdataset(f, chunks={"time": 1})[var]
        except Exception:
            da = open_mfdataset(f, chunks={})[var]
        if "vertical" in da.cf:
            da = check_lev(da)
        dsets.append(da)
    dsets += list(get_vc(ref_ds))
    return xr.merge(dsets, compat="override", join="override")


def get_gfile(scratch=None, **kwargs):
    if "df" in kwargs:
        df = kwargs["df"]
    else:
        # create a file catalog containing files and their
        # start and end times.
        files = create_catalog(**kwargs)
        data = dask.compute(files)
        df = create_df(data[0])
    return GFile(df=df, scratch=scratch)
