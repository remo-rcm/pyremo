import glob
import os
from datetime import timedelta as td
from pathlib import Path
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
    return {f: get_min_max_time(f) for f in files}


def get_files(dirs):
    """Get files from a list of directories"""
    if not isinstance(dirs, list):
        dirs = [dirs]
    files = []
    for d in dirs:
        d = Path(d).resolve()
        if d.is_dir():
            files += glob.glob(str(d / "*.nc"))
        else:
            files += glob.glob(str(d))
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
    """Scan directories for files and add time min and time max"""
    files = {v: get_files(d) for v, d in args.items()}
    for v, fs in files.items():
        files[v] = get_times_from_files(fs)
    return files


def to_cfdatetime(time, calendar="standard"):
    if isinstance(time, pd.Timestamp):
        # convert to string to handle precision issues with microseconds
        time = time.strftime("%Y-%m-%d %H:%M:%S")
    try:
        return xr.cftime_range(start=time, end=time, calendar=calendar)[0]
    except Exception:
        return np.nan


def search_df(df, **kwargs):
    """
    Search dataframe by arbitrary conditions.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to search.
    **kwargs : dict
        Arbitrary conditions to filter the dataframe. The keys are column names and the values are the conditions.

    Returns
    -------
    pandas.DataFrame
        The filtered dataframe.
    """
    condition_list = [
        (
            f"(df['{key}'].isin({repr(item)}))"
            if isinstance(item, list)
            else f"(df['{key}'] == {repr(item)})"
        )
        for key, item in kwargs.items()
    ]
    conditions = " & ".join(condition_list)
    return df[eval(conditions)]


def get_var_by_time(df, datetime=None, **kwargs):
    """
    Get variables by time.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to search.
    datetime : datetime, optional
        The datetime to filter the dataframe.
    **kwargs : dict
        Arbitrary conditions to filter the dataframe. The keys are column names and the values are the conditions.

    Returns
    -------
    pandas.DataFrame
        The filtered dataframe.
    """
    df = search_df(df, **kwargs)
    if datetime:
        df = df[(datetime >= df.time_min) & (datetime <= df.time_max)]
    return df


def get_sst_times(date):
    """
    Get daily dates from which the SST is interpolated in time.

    Parameters
    ----------
    date : datetime
        The date to get the SST times.

    Returns
    -------
    list of cftime.datetime
        The list of SST times.
    """
    if date.hour == 12:
        return (date,)
    days = (date + td(days=i) for i in range(-2, 3))
    return [
        cfdt.datetime(d.year, d.month, d.day, 12, calendar="standard") for d in days
    ]


def cdo_call(self, options="", op="", input="", output="temp", print_command=True):
    """
    Call CDO (Climate Data Operators) with the given options and operation.

    Parameters
    ----------
    options : str, optional
        The options to pass to CDO.
    op : str, optional
        The CDO operation to perform.
    input : str, optional
        The input file for the CDO operation.
    output : str, optional
        The output file for the CDO operation.
    print_command : bool, optional
        Whether to print the CDO command.

    Returns
    -------
    None
    """
    pass


class CFModelSelector:
    """
    CF Model Selector

    The CF model selector class simply selects files by a certain timestep.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe containing the file information.
    calendar : str, optional
        The calendar type.
    **kwargs : dict
        Arbitrary conditions to filter the dataframe.

    Attributes
    ----------
    calendar : str
        The calendar type.
    tempfiles : list
        The list of temporary files.
    df : pandas.DataFrame
        The filtered dataframe.
    """

    def __init__(self, df, calendar="standard", **kwargs):
        df = df.copy()
        if kwargs:
            df = search_df(df, **kwargs)
        self.calendar = calendar
        self.tempfiles = []
        self.df = self._update_time(df)

    def _update_time(self, df):
        """
        Update the time columns in the dataframe.

        Parameters
        ----------
        df : pandas.DataFrame
            The dataframe to update.

        Returns
        -------
        pandas.DataFrame
            The updated dataframe.
        """
        df["time_min"] = df["time_min"].apply(to_cfdatetime, calendar=self.calendar)
        df["time_max"] = df["time_max"].apply(to_cfdatetime, calendar=self.calendar)
        return df

    def get_file(self, datetime=None, **kwargs):
        """
        Get the file for the given datetime and conditions.

        Parameters
        ----------
        datetime : datetime, optional
            The datetime to filter the dataframe.
        **kwargs : dict
            Arbitrary conditions to filter the dataframe. The keys are column names and the values are the conditions.

        Returns
        -------
        str
            The path to the file.

        Raises
        ------
        FileNotFoundError
            If no file is found for the given conditions.
        """
        if datetime:
            datetime = to_cfdatetime(datetime, self.calendar)
        sel = get_var_by_time(self.df, datetime=datetime, **kwargs)
        if len(sel.index) > 1:
            return list(sel.path)
        if sel.empty:
            raise FileNotFoundError(f"no file found: {kwargs}, datetime: {datetime}")
        return sel.iloc[0].path


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
    dsets += list(get_vc(ref_ds))
    output = xr.merge(dsets, compat="override", join="override")
    output.attrs = ref_ds.attrs
    return output


def get_gcm_dataset(
    ds,
    ref_ds=None,
    tos=None,
    attrs=None,
    use_cftime=True,
    invertlev=None,
    time_range=None,
):
    """Creates a global dataset ready for preprocessing.

    This function creates a homogenized global dataset. If necessary,
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
        ds = open_datasets(ds, ref_ds, time_range=time_range)
    else:
        ds = ds.copy()
        ds["akgm"], ds["bkgm"] = get_vc(ds, invertlev)
        ds = check_lev(ds, invertlev)

    if tos is not None:
        print("adding sst...")
        ds["tos"] = map_sst(tos, ds)
    # ensure correct units
    ds = convert_units(ds)
    if "sftlf" in ds:
        ds["sftlf"] = np.around(ds.sftlf)
    if attrs is None:
        attrs = ds.attrs
    if use_cftime is True:
        return ds.assign_coords(
            time=ds.time.convert_calendar(ds.time.dt.calendar, use_cftime=True).time
        )
    return ds


def gfile(*args, **kwargs):
    warn(
        "The 'gfile' function is deprecated and will be removed in a future version. "
        "Please use 'get_gcm_dataset' instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return get_gcm_dataset(*args, **kwargs)


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
                self.scratch = Path(os.environ["SCRATCH"]) / ".cf-selector"
            except Exception:
                pass

        self.tos_regridder = None

    def get_files(self, variables, datetime=None, **kwargs):
        if datetime:
            datetime = to_cfdatetime(datetime, self.calendar)
        return {
            var: self.selector.get_file(variable_id=var, datetime=datetime, **kwargs)
            for var in variables
        }

    def extract_timestep(self, datetime=None, **kwargs):
        pass

    def extract_dynamic_timesteps(self, datetime=None, **kwargs):
        """Extract timesteps for a certain date from input files."""
        datetime = to_cfdatetime(datetime, self.calendar)
        files = self.get_files(self.dynamics, datetime=datetime, **kwargs)
        if datetime is None:
            return files
        cdo = Cdo(tempdir=self.scratch) if self.scratch else Cdo()
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
        files = {}
        for t in times:
            try:
                f = self.selector.get_file(variable_id=self.sst, datetime=t)
                files[t] = f
            except FileNotFoundError:
                warn(f"sst not found for {datetime}, will extrapolate...")
                pass

        cdo = Cdo(tempdir=self.scratch) if self.scratch else Cdo()
        sst_extract = [
            cdo.seldate(
                t.strftime(cdo_datetime_format), input=f, returnXDataset=True
            ).load()
            for t, f in files.items()
        ]
        sst_da = xr.merge(sst_extract)[self.sst]
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
        """Regrid tos to atmospheric grid."""
        import xesmf as xe

        attrs = da.attrs
        atmo_grid = atmo_grid.copy()
        atmo_grid["mask"] = ~(atmo_grid.sftlf > 0).squeeze(drop=True)

        ds = da.to_dataset()

        if not self.tos_regridder:
            ds["mask"] = ~ds.tos.isnull().squeeze(drop=True)
            self.tos_regridder = xe.Regridder(
                ds, atmo_grid, method="nearest_s2d", extrap_method="nearest_s2d"
            )

        out = self.tos_regridder(da)
        out.attrs = attrs
        return out

    def gfile(self, datetime, sst=True, **kwargs):
        """Creates a gfile from CF input data."""
        gds = self.extract_data(datetime=datetime, **kwargs)
        if sst:
            gds[self.sst] = self.get_sst(datetime, gds)
        if "variable_id" in gds.attrs:
            del gds.attrs["variable_id"]
        return gfile(gds)


def get_gcm_gfile(scratch=None, **kwargs):
    """
    Create a GFile object for GCM data.

    This function creates a GFile object which is used to handle and process
    global climate model (GCM) data. It can either take a pre-existing dataframe
    or create one by scanning directories for files and extracting time information.

    Parameters
    ----------
    scratch : str or None, optional
        Path to a scratch directory for temporary files. If None, the default
        scratch directory is used.
    **kwargs : dict
        This should be a dictionary of cf variable names as keys and directories
        of netCDF files as values. These are scanned to be used for the preprocessing.

    Returns
    -------
    GFile
        An instance of the GFile class containing the processed GCM data.
    """
    if "df" in kwargs:
        df = kwargs["df"]
    else:
        files = create_catalog(**kwargs)
        data = dask.compute(files)
        df = create_df(data[0])
    return GFile(df=df, scratch=scratch)


def get_gfile(*args, **kwargs):
    warn(
        "The 'get_gfile' function is deprecated and will be removed in a future version. "
        "Please use 'get_gcm_gfile' instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return get_gcm_gfile(*args, **kwargs)
