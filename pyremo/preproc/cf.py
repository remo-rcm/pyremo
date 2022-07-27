import os

# import subprocess
# import tempfile
from datetime import timedelta as td

import cftime as cfdt
import numpy as np
import pandas as pd
import xarray as xr
from cdo import Cdo

from .constants import lev_i
from .core import check_lev, convert_units, get_vc2, horizontal_dims, open_mfdataset

cdo_exe = "cdo"
default_catalog = "/work/ik1017/Catalogs/dkrz_cmip6_disk.csv.gz"

cdo_datetime_format = "%Y-%m-%dT%H:%M:%S"


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


class CFModelSelector:
    def __init__(self, df=None, scratch=None, calendar="standard", **kwargs):
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
            raise Exception("no file found")
        return sel.iloc[0].path

    # own cdo call here due to https://github.com/Try2Code/cdo-bindings/issues/34
    # def _cdo_call(self, options="", op="", input="", output="temp", print_command=True):
    #     if output is None:
    #         output = ""
    #     elif output == "temp":
    #         output = tempfile.TemporaryDirectory(dir=self.scratch).name
    #         self.tempfiles.append(output)
    #     if isinstance(input, list):
    #         input = " ".join(input)
    #     call = "{} {} {} {} {}".format(cdo_exe, options, op, input, output)
    #     if print_command is True:
    #         print(call)
    #     stdout = subprocess.Popen(
    #         call, shell=True, stdout=subprocess.PIPE
    #     ).stdout.read()
    #     if output:
    #         return output
    #     return stdout
    def _cdo_call(self, options="", op="", input="", output="temp", print_command=True):
        cdo = Cdo(tempdir=self.scratch)
        getattr(cdo, op)(options=options, input=input)


def gfile(ds, ref_ds=None, tos=None, time_range=None, attrs=None):
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
    #    if tos is not None:
    #        ds["tos"] = map_sst(tos, ds.sel(time=time_range))
    ds = ds.rename({"lev": lev_i})
    ds = convert_units(ds)
    if "sftlf" in ds:
        ds["sftlf"] = np.around(ds.sftlf)
    if attrs is None:
        attrs = ds.attrs
    return ds


class GFile:

    dynamics = ["ta", "ua", "va", "ps", "hus"]
    fx = ["orog", "sftlf"]
    sst = "tos"
    all_vars = dynamics + fx  # + sst

    def __init__(self, df=None, scratch=None, **kwargs):
        self.selector = CFModelSelector(df=df, **kwargs)
        if scratch is None:
            try:
                scratch = os.path.join(os.environ["SCRATCH"], ".cf-selector")
            except Exception:
                pass
        self.cdo = Cdo(tempdir=scratch, logging=True)
        self.regridder = None

    def get_files(self, variables, datetime, **kwargs):
        datetime = to_cfdatetime(datetime)
        files = {
            var: self.selector.get_file(variable_id=var, datetime=datetime, **kwargs)
            for var in variables
        }
        return files

    def extract_timestep(self, datetime=None, **kwargs):
        pass

    def extract_dynamic_timesteps(self, datetime=None, **kwargs):
        files = self.get_files(self.dynamics, datetime=datetime, **kwargs)
        if datetime is None:
            return files
        return {v: self.cdo.seldate(datetime, input=f) for v, f in files.items()}

    def extract_files(self, datetime, **kwargs):
        files = self.extract_dynamic_timesteps(datetime, **kwargs)
        files.update(self.get_files(self.fx, datetime=None, **kwargs))
        # files.update(self.get_files(self.sst, datetime=datetime, **kwargs))
        return files

    def get_sst(self, datetime, atmo_grid):
        datetime = to_cfdatetime(datetime)
        times = get_sst_times(datetime)
        files = {
            t: self.selector.get_file(variable_id=self.sst, datetime=t) for t in times
        }
        sst_extract = [
            self.cdo.seldate(t.strftime(cdo_datetime_format), input=f)
            for t, f in files.items()
        ]
        sst_ds = xr.open_mfdataset(sst_extract, use_cftime=True)
        sst_ds = sst_ds[self.sst].interp(time=datetime.strftime(cdo_datetime_format))
        return self.regrid_to_atmosphere(sst_ds, atmo_grid)

    def regrid_to_atmosphere(self, ds, atmo_grid):
        import xesmf as xe

        atmo_grid = atmo_grid.copy()
        if self.regridder is None:
            print("creating regridder")
            self.regridder = xe.Regridder(ds, atmo_grid, method="nearest_s2d")
        ds = self.regridder(ds)
        return ds

    def get_gfile(self, datetime, sst=True, **kwargs):
        files = self.extract_files(datetime=datetime, **kwargs)
        gds = xr.open_mfdataset(files.values(), compat="override", join="override")
        if sst is True:
            sst_ds = self.get_sst(datetime, gds)
            sst_ds.name = self.sst
            gds = gds.merge(sst_ds)
        return gds


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
