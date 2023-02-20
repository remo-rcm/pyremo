import glob
import os

# import subprocess
# import tempfile
from datetime import timedelta as td
from os import path as op

import cftime as cfdt
import dask
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from cdo import Cdo

# from .constants import lev_i
from .core import check_lev, convert_units, get_vc2, horizontal_dims, open_mfdataset

cdo_exe = "cdo"
default_catalog = "/work/ik1017/Catalogs/dkrz_cmip6_disk.csv.gz"

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


# @dask.delayed
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
    # time = str(int(time))
    try:
        return xr.cftime_range(start=time, end=time, calendar=calendar)[0]
    except Exception:
        return np.nan


def search_df(df, **kwargs):
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
    if datetime is not None and len(df) > 1:
        df = df[(datetime >= df.time_min) & (datetime <= df.time_max)]
    return df


def cdo_call(self, options="", op="", input="", output="temp", print_command=True):
    pass


class CFModelSelector:
    def __init__(self, df=None, calendar="standard", **kwargs):
        if df is None:
            df = pd.read_csv(default_catalog)
        df = df.copy()
        if kwargs:
            df = search_df(df, **kwargs)
        self.calendar = calendar
        self.tempfiles = []
        self.df = self._update_time(df)

    def __repr__(self):
        return repr(self._group())

    def _repr_html_(self):
        return self._group()._repr_html_()

    def _group(self):
        groups = ["source_id", "member_id", "experiment_id", "table_id"]
        return self.df.groupby(groups)[
            [
                "variable_id",
                #  "member_id",
                "institution_id",
                #  "table_id",
                "activity_id",
            ]
        ].agg(["unique"])

    def _update_time(self, df):
        df["time_min"] = df["time_min"].apply(to_cfdatetime, calendar=self.calendar)
        df["time_max"] = df["time_max"].apply(to_cfdatetime, calendar=self.calendar)
        return df

    def get_file(self, datetime=None, **kwargs):
        sel = get_var_by_time(self.df, datetime=datetime, **kwargs)
        if len(sel.index) > 1:
            return list(sel.path)
            # raise Exception("file selection is not unique")
        if sel.empty:
            raise Exception("no file found: {}, date: {}".format(kwargs, datetime))
        return sel.iloc[0].path

    # def _cdo_call(self, options="", op="", input="", output="temp", print_command=True):
    #    cdo = Cdo(tempdir=self.scratch)
    #    getattr(cdo, op)(options=options, input=input)


def gfile(ds, ref_ds=None, tos=None, time_range=None, attrs=None, use_cftime=True):
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

    tos : xarray.Dataset
        Sea surface dataset.

    time_rage :
        The common time range from the input and sst that should be used.

    attrs:
        Global attributes for the output dataset. If ``attrs=None``, the global
        attributes from ``ref_ds`` are used.

    Returns
    -------
    gfile : xarray.Dataset
        Global dataset ready for preprocessing.

    """
    if isinstance(ds, dict):
        ds = open_datasets(ds, ref_ds, time_range)
        if time_range is None:
            time_range = ds.time
    else:
        ds = ds.copy()
        if time_range is None:
            time_range = ds.time
        ds = ds.sel(time=time_range)
        ds["akgm"], ds["bkgm"] = get_vc2(ds)
        ds = check_lev(ds)
    # if tos is not None:
    #    ds["tos"] = map_sst(tos, ds.sel(time=time_range))
    # ds = ds.rename({"lev": lev_i})
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
        # self.cdo = Cdo(tempdir=self.scratch)
        self.regridder = None

    def get_files(self, variables, datetime, **kwargs):
        datetime = to_cfdatetime(datetime, self.calendar)
        files = {
            var: self.selector.get_file(variable_id=var, datetime=datetime, **kwargs)
            for var in variables
        }
        return files

    def extract_timestep(self, datetime=None, **kwargs):
        pass

    def extract_dynamic_timesteps(self, datetime=None, **kwargs):
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
        # files.update(self.get_files(self.sst, datetime=datetime, **kwargs))
        return xr.merge(files.values(), compat="override", join="override")
        # return files

    def get_sst(self, datetime, atmo_grid=None):
        datetime = to_cfdatetime(datetime, self.calendar)
        times = get_sst_times(datetime)
        files = {
            t: self.selector.get_file(variable_id=self.sst, datetime=t) for t in times
        }
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
            sst_da = sst_da.interp(time=datetime.strftime(cdo_datetime_format))
        if atmo_grid is None:
            return sst_da
        return self.regrid_to_atmosphere(sst_da.squeeze(drop=True), atmo_grid)

    def regrid_to_atmosphere(self, da, atmo_grid):
        import xesmf as xe

        attrs = da.attrs
        atmo_grid = atmo_grid.copy()
        atmo_grid["mask"] = ~(atmo_grid.sftlf > 0).squeeze(drop=True)
        # if self.regridder is None:
        ds = da.to_dataset()
        ds["mask"] = ~ds.tos.isnull().squeeze(drop=True)
        print("creating regridder")
        self.regridder = xe.Regridder(
            ds, atmo_grid, method="nearest_s2d", extrap_method="nearest_s2d"
        )
        da = self.regridder(ds.tos)
        da.attrs = attrs
        return da

    def gfile(self, datetime, sst=True, **kwargs):
        gds = self.extract_data(datetime=datetime, **kwargs)
        if sst is True:
            gds[self.sst] = self.get_sst(datetime, gds)
        if "variable_id" in gds.attrs:
            del gds.attrs["variable_id"]
        return gfile(gds)
        # gds = convert_units(gds)
        # if "sftlf" in gds:
        #    gds["sftlf"] = np.around(gds.sftlf)
        # return gds


def get_sst_times(dt):
    cal = dt.calendar
    if dt.hour == 12:
        return (dt,)
    if dt.hour > 12:
        dt1 = dt + td(days=1)
    else:
        dt1 = dt
        dt = dt + td(days=-1)
    return (
        cfdt.datetime(dt.year, dt.month, dt.day, 12, calendar=cal),
        cfdt.datetime(dt1.year, dt1.month, dt1.day, 12, calendar=cal),
    )


def open_datasets(datasets, ref_ds=None, time_range=None):
    """Creates a virtual gfile"""
    if ref_ds is None:
        try:
            ref_ds = open_mfdataset(datasets["ta"])
        except Exception:
            raise Exception("ta is required in the datasets dict if no ref_ds is given")
    lon, lat = horizontal_dims(ref_ds)
    # ak_bnds, bk_bnds = get_ab_bnds(ref_ds)
    if time_range is None:
        time_range = ref_ds.time
    dsets = []
    for var, f in datasets.items():
        try:
            da = open_mfdataset(f, chunks={"time": 1})[var]
            da = da.sel(time=time_range)
        except Exception:
            da = open_mfdataset(f, chunks={})[var]
        if "vertical" in da.cf:
            da = check_lev(da)
        dsets.append(da)
    dsets += list(get_vc2(ref_ds))
    return xr.merge(dsets, compat="override", join="override")


def get_gfile(scratch=None, **kwargs):
    if "df" in kwargs:
        df = kwargs["df"]
    else:
        files = create_catalog(**kwargs)
        data = dask.compute(files)
        df = create_df(data[0])
    return GFile(df=df, scratch=scratch)
