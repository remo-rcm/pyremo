import os
from os import path as op

import cordex as cx
import numpy as np
import xarray as xr

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
    for dim in da.dims:
        if "lon" in dim:
            lon_dim = dim
        if "lat" in dim:
            lat_dim = dim
    return (lon_dim, lat_dim)


def get_grid(domain_info):
    return cx.create_dataset(**domain_info)


def encoding(da, missval=None):
    result = {}
    result.update(encode_missval(da, missval))
    return result


def encode_missval(da, missval=None):
    if missval is None:
        missval = 1.0e20
    if np.isnan(da.values).any():
        return {"_FillValue": missval}
    else:
        return {"_FillValue": None}


def update_attrs(ds):
    """Update attributes of all variables for CF compliance"""
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
    if template is None:
        template = "x{}x{}.nc"
    return template.format(expid, date.strftime("%Y%m%d%H"))


def write_forcing_file(
    ds,
    path=None,
    expid="000000",
    template=None,
    tempfiles=None,
    missval=1.0e20,
    **kwargs
):
    if path is None:
        path = "./"
    if not os.path.isdir(path):
        os.makedirs(path)
    if template is None:
        template = "a{}a{}.nc"
    fname = op.join(path, get_filename(ds.time.data[0], expid, template))
    for v in ds.data_vars:
        var = ds[v]
        var.encoding = encoding(var, missval)
        # expand time dim is neccessary for REMO input forcing
        if v not in static and "time" not in var.dims:
            ds[v] = var.expand_dims("time")
    ds.to_netcdf(fname, **kwargs)
    return fname


def to_netcdf(
    ads,
    path="",
    expid="000000",
    template=None,
    tempfiles=None,
    missval=1.0e20,
    **kwargs
):
    """write dataset to netcdf

    by default, each timestep goes into a separate output file

    """
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
    for var, da in ds.items():
        if var in expand_time:
            ds[var] = da.expand_dims("time")
        if var in fillvars:
            ds[var].encoding["_FillValue"] = missval
        else:
            ds[var].encoding["_FillValue"] = None
