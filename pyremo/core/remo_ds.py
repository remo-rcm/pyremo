# -*- coding: utf-8 -*-
# flake8: noqa
"""Remo dataset module.

This module contains functions to work with REMO datasets.

"""

# flake8: noqa
from . import cal, codes
from .cal import parse_dates


def preprocess(ds, use_cftime=False):
    """preprocessing for opening with xr.open_mfdataset

    This function can be used as the preprocess function for
    opening a REMO dataset with xr.open_mfdataset. The function
    will update meta information according to the REMO code table
    and also parse the time axis if it contains absolute times.

    """
    ds = update_meta_info(ds)
    try:
        return parse_dates(ds, use_cftime)
    except:
        return ds


def open_remo_mfdataset(filenames, update_meta=False, parse_dates=False):
    import xarray as xr

    ds = xr.open_mfdataset(filenames)
    if update_meta:
        ds = update_meta_infos(ds)
    if parse_dates:
        ds = parse_dates(ds)
    return ds


def open_remo_dataset(
    filename,
    options="-f nc4",
    update_meta=False,
    returnX=True,
    inplace=False,
    parse_time=False,
    **kwargs
):
    """Read a REMO dataset.

    Read in a REMO dataset into xarray.Dataset or netCDF4.Dataset from IEG or NetCDF and
    provide meta data.

    Parameters
    ----------
    filename : str
        Filename of IEG or NetCDF file.
    update_meta: bool
        Update variable meta information of the dataset from REMO tables.
    returnX : bool
        Return an xarray.Dataset. If False, use netCDF4 Dataset.
    inplace: bool
        Update meta info on disk, only useful for netCDF4 Datasets.
    parse_time: bool
        Parse absolute time axis into datetime objects.

    Returns
    -------
    dataset
        Returns an xarray.Dataset or netCDF4.Dataset.

    """
    format = _get_fileformat(filename)
    # print(format)
    if "NetCDF" in format and not options:
        ds = _read_nc_dataset(filename, returnX=returnX, inplace=inplace, **kwargs)
    elif "IEG" in format:
        if returnX:
            ds = _read_with_cdo(filename, options, returnXDataset=True)
        else:
            ds = _read_with_cdo(filename, options, returnCdf=not inplace)
            if inplace:
                ds = _read_nc_dataset(ds, returnX=False, inplace=inplace)
    else:
        ds = _read_nc_dataset(filename, returnX, **kwargs)
    if update_meta:
        ds = update_meta_info(ds)
    if parse_time is True:
        ds = cal.parse_dates(ds)
    return ds


def _read_nc_dataset(filename, returnX=True, inplace=False, **kwargs):
    """Use xarray or netCDF4 to read NetCDF."""
    if returnX:
        import xarray as xr

        if type(filename) is list:
            return xr.open_mfdataset(filename, **kwargs)
        return xr.open_dataset(filename, **kwargs)
    else:
        from netCDF4 import Dataset, MFDataset

        if inplace:
            mode = "a"
        else:
            mode = "r"
        if type(filename) is list:
            return MFDataset(filename, mode="r")
        return Dataset(filename, mode="a")


def _read_with_cdo(filename, options="", **kwargs):
    """uses cdo to read unknown file format."""
    from cdo import Cdo

    return Cdo().copy(options=options, input=filename, **kwargs)


def _get_fileformat(filename):
    """"""
    try:
        from cdo import Cdo

        return Cdo().showformat(input=filename)[0]
    except:
        return "Unknown"


def update_meta_info(ds, id=None):
    """Updates meta info of a dataset.

    Updates variable names and attributes in an xarray Dataset
    based on Remo table.

    Parameters
    ----------
    ds : dataset
        The dataset in which meta info should be updated.

    Returns
    -------
    dataset
        The dataset with updated meta information.

    """
    meta = {var: _get_meta_info(ds[var]) for var in ds.variables}
    renames = {var: info["variable"] for var, info in meta.items() if info if not None}
    ds = _update_attrs(ds, meta)
    ds = _rename_ds(ds, renames)
    return ds


def _get_meta_info(da, id=None):
    """Get meta info for a netcdf variable."""
    if id:
        attrs = codes.get_dict(id)
    else:
        try:
            attrs = codes.get_dict(da.name)
        except:
            try:
                attrs = codes.get_dict(da.code)
            except:
                attrs = None
    return attrs


def _rename_ds(ds, renames):
    """Rename variables in a dataset."""
    try:  # xarray
        return ds.rename(renames)
    except:  # netCDF4
        for old, new in renames.items():
            if old != new:
                ds.renameVariable(old, new)
        return ds


def _update_attrs(ds, meta):
    """Update variable attributes in a dataset."""
    for var, info in meta.items():
        if info:
            filter_info = {key: value for key, value in info.items() if value}
            # print(filter_info)
            try:  # xarray
                ds[var].attrs.update(filter_info)
            except:  # netCDF4
                ds[var].setncatts(filter_info)
    return ds
