import os
from os import path as op
from warnings import warn
import cordex as cx
import numpy as np
import xarray as xr
import pandas as pd

import pyremo as pr

# variables that should have a mask with fill values
fillvars = ["TSW", "SEAICE", "TSI"]
vcs = [
    "hyai",
    "hybi",
    "hyam",
    "hybm",
    "akgm",
    "bkgm",
    "ak",
    "bk",
    "rotated_latitude_longitude",
]
static = vcs + ["rotated_latitude_longitude"]


def horizontal_dims(da):
    """
    Identify horizontal dimensions in the dataset.

    Parameters
    ----------
    da : xarray.DataArray
        Input data array.

    Returns
    -------
    tuple
        Tuple containing longitude and latitude dimensions.
    """
    for dim in da.dims:
        if "lon" in dim:
            lon_dim = dim
        if "lat" in dim:
            lat_dim = dim
    return (lon_dim, lat_dim)


def get_grid(domain_info):
    """
    Create a dataset based on domain information.

    Parameters
    ----------
    domain_info : dict
        Dictionary containing domain information.

    Returns
    -------
    xarray.Dataset
        Created dataset.
    """
    return cx.create_dataset(**domain_info)


def encoding(da, missval=None):
    """
    Generate encoding for a data array.

    Parameters
    ----------
    da : xarray.DataArray
        Input data array.
    missval : float, optional
        Missing value, by default None.

    Returns
    -------
    dict
        Encoding dictionary.
    """
    result = {}
    result.update(encode_missval(da, missval))
    return result


def encode_missval(da, missval=None):
    """
    Encode missing values in a data array.

    Parameters
    ----------
    da : xarray.DataArray
        Input data array.
    missval : float, optional
        Missing value, by default None.

    Returns
    -------
    dict
        Dictionary with '_FillValue' key.
    """
    if missval is None:
        missval = 1.0e20
    if np.isnan(da.values).any():
        return {"_FillValue": missval}
    else:
        return {"_FillValue": None}


def update_attrs(ds):
    """
    Update attributes of all variables for CF compliance.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset.

    Returns
    -------
    xarray.Dataset
        Dataset with updated attributes.
    """
    for var, da in ds.items():
        try:
            attrs = pr.codes.get_dict(var)
            da.attrs = {}
            for attr, value in attrs.items():
                if value is not None:
                    da.attrs[attr] = value

            da.attrs["grid_mapping"] = "rotated_latitude_longitude"
            da.attrs["coordinates"] = "lon lat"

        except Exception:
            pass
    return ds


def get_filename(date, expid="000000", template=None):
    """
    Generate a filename based on date and experiment ID.

    Parameters
    ----------
    date : datetime
        Date for the filename.
    expid : str, optional
        Experiment ID, by default "000000".
    template : str, optional
        Filename template, by default None.

    Returns
    -------
    str
        Generated filename.
    """
    if template is None:
        template = "a{expid}a{date:%Y%m%d%H}.nc"
    return template.format(expid=expid, date=date)


def write_forcing_file(
    ds, path=None, expid="000000", template=None, missval=1.0e20, **kwargs
):
    """
    Write a forcing file.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to write.
    path : str, optional
        Output path, by default None.
    expid : str, optional
        Experiment ID, by default "000000".
    template : str, optional
        Filename template, by default None.
    missval : float, optional
        Missing value, by default 1.0e20.

    Returns
    -------
    str
        Path to the written file.
    """
    if path is None:
        path = "./"
    fname = op.join(path, get_filename(ds.time.item(), expid, template))
    for v in ds.data_vars:
        var = ds[v]
        if v in fillvars:
            var.encoding = {"_FillValue": missval}
        else:
            var.encoding = {"_FillValue": None}
        # expand time dim is necessary for REMO input forcing
        if v not in static and "time" not in var.dims:
            ds[v] = var.expand_dims("time")
    ds.to_netcdf(fname, **kwargs)
    ds.close()
    return fname


def to_netcdf(
    ads,
    path="",
    expid="000000",
    template=None,
    tempfiles=None,
    missval=1.0e20,
    **kwargs,
):
    """
    Write dataset to netCDF.

    By default, each timestep goes into a separate output file.

    Parameters
    ----------
    ads : xarray.Dataset
        Input dataset.
    path : str, optional
        Output path, by default "".
    expid : str, optional
        Experiment ID, by default "000000".
    template : str, optional
        Filename template, by default None.
    tempfiles : list, optional
        List of temporary files, by default None.
    missval : float, optional
        Missing value, by default 1.0e20.

    Returns
    -------
    list
        List of paths to the written files.
    """
    warn(
        "The 'to_netcdf' function is deprecated and will be removed in a future version. "
        "Please use the 'write_forcing_file' function instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    if not os.path.isdir(path):
        os.makedirs(path)
    expand_time = [var for var, da in ads.items() if "time" in da.dims]
    if template is None:
        template = "a{}a{}.nc"
    dates, datasets = zip(*ads.groupby("time"))
    paths = [os.path.join(path, get_filename(date, expid, template)) for date in dates]
    dsets = []
    # expand time dimension only for variables not coordinates.
    for ds in datasets:
        for var, da in ds.items():
            if var in expand_time:
                ds[var] = da.expand_dims("time")
            if var in fillvars:
                ds[var].encoding["_FillValue"] = missval
            else:
                ds[var].encoding["_FillValue"] = None
        dsets.append(ds)
    # dsets = [dset.expand_dims('time') for dset in datasets]
    xr.save_mfdataset(dsets, paths, **kwargs)
    if tempfiles is not None:
        for f in tempfiles:
            os.remove(f)
    return paths


def to_tar(files, tar_file, mode="w"):
    """
    Create a tar file from a list of files.

    Parameters
    ----------
    files : list
        List of files to include in the tar file.
    tar_file : str
        Path to the output tar file.
    mode : str, optional
        Mode for opening the tar file, by default "w".

    Returns
    -------
    str
        Path to the created tar file.
    """
    import tarfile

    try:
        from tqdm import tqdm
    except Exception:

        def tqdm(x):
            return x

    tf = tarfile.open(tar_file, mode=mode)
    for f in tqdm(files, desc="creating tarfile"):
        tf.add(f, arcname=os.path.basename(f), recursive=False)
    tf.close()
    return tar_file


def encode(ds, expand_time, missval=1.0e20):
    """
    Encode dataset with missing values and expand time dimension.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset.
    expand_time : list
        List of variables to expand time dimension.
    missval : float, optional
        Missing value, by default 1.0e20.

    Returns
    -------
    xarray.Dataset
        Encoded dataset.
    """
    for var, da in ds.items():
        if var in expand_time:
            ds[var] = da.expand_dims("time")
        if var in fillvars:
            ds[var].encoding["_FillValue"] = missval
        else:
            ds[var].encoding["_FillValue"] = None
    return ds


def datelist(startdate, enddate, freq, inclusive="both", **kwargs):
    return pd.date_range(startdate, enddate, freq=freq, inclusive=inclusive, **kwargs)


def ensure_dir(path):
    if not os.path.exists(path):
        # logger.debug(f"Creating directory {path}")
        try:
            os.makedirs(path)
        except Exception as e:
            warn(e)
