"""core module for preprocessing

This module wraps the pyintorg interfaces into xr.apply_ufunc.

"""

import warnings

import cf_xarray as cfxr
import numpy as np
import xarray as xr

from .constants import const, lev, lev_i
from .utils import horizontal_dims

xr.set_options(keep_attrs=True)


def open_mfdataset(
    files,
    use_cftime=True,
    parallel=True,
    data_vars="minimal",
    chunks={},
    coords="minimal",
    compat="override",
    drop=None,
    **kwargs,
):
    """optimized function for opening CMIP6 6hrLev 3d datasets

    based on https://github.com/pydata/xarray/issues/1385#issuecomment-561920115

    """

    def drop_all_coords(ds):
        # ds = ds.drop(drop)
        return ds.reset_coords(drop=True)

    ds = xr.open_mfdataset(
        files,
        parallel=parallel,
        decode_times=False,
        combine="by_coords",
        preprocess=drop_all_coords,
        decode_cf=False,
        chunks=chunks,
        data_vars=data_vars,
        coords="minimal",
        compat="override",
        **kwargs,
    )
    return xr.decode_cf(ds, use_cftime=use_cftime)


def get_akbkem(vc):
    """create vertical coordinate dataset"""
    akbk = vc.to_xarray().drop("index")
    # bkem = pr.tables.vc.tables['vc_27lev']
    akem = akbk.ak.swap_dims({"index": lev_i})
    bkem = akbk.bk.swap_dims({"index": lev_i})
    akhem = (0.5 * (akbk.ak[:-1] + akbk.ak[1:])).swap_dims({"index": lev})
    bkhem = (0.5 * (akbk.bk[:-1] + akbk.bk[1:])).swap_dims({"index": lev})
    akem[lev_i] = xr.DataArray(np.arange(1, akem.size + 1), dims=lev_i)
    bkem[lev_i] = xr.DataArray(np.arange(1, bkem.size + 1), dims=lev_i)
    akhem[lev] = xr.DataArray(np.arange(1, akhem.size + 1), dims=lev, name="akh")
    bkhem[lev] = xr.DataArray(np.arange(1, bkhem.size + 1), dims=lev, name="bkh")
    akhem.name = "akh"
    bkhem.name = "bkh"
    return xr.merge([akem, bkem, akhem, bkhem])


def get_ab_bnds(ds):
    ak_valid = ["ap_bnds", "a_bnds", "hyai"]
    bk_valid = ["b_bnds", "hybi"]
    ak_bnds = None
    bk_bnds = None
    for ak_name in ak_valid:
        if ak_name in ds:
            ak_bnds = ds[ak_name]
            print("using {} for akgm".format(ak_name))
    for bk_name in bk_valid:
        if bk_name in ds:
            bk_bnds = ds[bk_name]
            print("using {} for bkgm".format(bk_name))
    return ak_bnds, bk_bnds


def get_vc(ds, invert=None):
    """Reads the vertical hybrid coordinate from a dataset."""
    if ds.cf["vertical"].attrs.get("positive") == "down" and invert is None:
        invert = True
    ak_bnds, bk_bnds = get_ab_bnds(ds)
    if ak_bnds.ndim > 1:
        ak = cfxr.bounds_to_vertices(ak_bnds, bounds_dim="bnds")
        bk = cfxr.bounds_to_vertices(bk_bnds, bounds_dim="bnds")
    else:
        ak = ak_bnds
        bk = bk_bnds
    if invert is True:
        print("inverting vertical coordinates")
        ak = np.flip(ak)
        bk = np.flip(bk)
    ak.name = "akgm"
    bk.name = "bkgm"
    return ak, bk


def map_sst(tos, ref_ds, resample="6H", regrid=True):
    from datetime import timedelta as td

    import xesmf as xe

    try:
        tos = tos.to_dataset()
    except Exception:
        pass
    # tos_res = tos
    attrs = tos.tos.attrs
    tos_times = (ref_ds.time.min() - td(days=1), ref_ds.time.max() + td(days=1))
    tos = tos.sel(time=slice(tos_times[0], tos_times[1]))
    # return tos_res
    # tos = tos.resample(time=resample).interpolate("linear").chunk({"time": 1})
    tos = tos.resample(time=resample).interpolate("linear")
    tos = tos.sel(time=ref_ds.time)

    if regrid:
        ref_ds["mask"] = ~(ref_ds.sftlf > 0)
        tos["mask"] = ~tos.tos.isel(time=0).isnull().squeeze(drop=True)
        regridder = xe.Regridder(tos, ref_ds, "nearest_s2d")
        tos = regridder(tos.tos)
    tos.attrs.update(attrs)

    return tos


def convert_units(ds):
    """convert units of input data for use in the preprocessor"""
    try:
        if ds.sftlf.units == "%":
            print("converting sftlf units to fractional")
            attrs = ds.sftlf.attrs
            ds["sftlf"] = ds.sftlf * 0.01
            attrs["units"] = 1
            ds.sftlf.attrs = attrs
    except Exception:
        warnings.warn("sftlf has no units attribute, must be fractional.")
    try:
        if ds.tos.units == "degC":
            print("converting tos units to K")
            attrs = ds.tos.attrs
            ds["tos"] = ds.tos + const.absolute_zero
            attrs["units"] = "K"
            ds.tos.attrs = attrs
    except Exception:
        warnings.warn("tos has no units attribute, must be Kelvin!")
    try:
        if ds.orog.units == "m":
            print("converting orography to geopotential")
            attrs = ds.orog.attrs
            ds["orog"] = ds.orog * const.grav_const
            attrs["units"] = "m2 s-2"
            ds.orog.attrs = attrs
    except Exception:
        warnings.warn("orog has no units attribute, must be m2 s-2")
    return ds


def check_lev(ds, invert=None):
    """Check for order of levels and invert if neccessary"""

    if ds.cf["vertical"].attrs.get("positive") == "down" and invert is None:
        invert = True
    # if "vertical" in ds.cf and auto is True:
    #    positive = ds.cf["vertical"].attrs.get("positive", None)
    # else:
    #    positive = None
    # if positive is None and auto is True:
    #    warnings.warn("could not determine positive attribute of vertical axis.")
    #    return ds
    if invert is True:
        kwargs = {ds.cf["vertical"].name: ds.cf["vertical"][::-1]}
        print("inverting vertical axis")
        return ds.reindex(**kwargs)
    return ds


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

    time_range :
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
        ds["akgm"], ds["bkgm"] = get_vc(ds)
        ds = check_lev(ds)
    if tos is not None:
        ds["tos"] = map_sst(tos, ds.sel(time=time_range))
    # ds = ds.rename({"lev": lev_gm})
    ds = convert_units(ds)
    # if "sftlf" in ds:
    #    ds["sftlf"] = np.around(ds.sftlf)
    if attrs is None:
        attrs = ds.attrs
    return ds
