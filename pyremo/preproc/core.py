"""core module for preprocessing

This module wraps the pyintorg interfaces into xr.apply_ufunc.

"""

import xarray as xr
import numpy as np
import warnings

try:
    from pyintorg import interface as intf
except:
    print(
        "could not find pyintorg, you need this for preprocessing. Please consider installing it from https://git.gerics.de/python/pyintorg.git"
    )


lev_i = "lev_i"
lev = "lev"
lev_gm = "lev_gm"


class const:
    """constants used for unit conversion"""

    grav_const = 9.806805923
    absolute_zero = 273.5


def open_mfdataset(
    files,
    use_cftime=True,
    parallel=True,
    data_vars="minimal",
    chunks={"time": 1},
    coords="minimal",
    compat="override",
    drop=None,
    **kwargs
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
        **kwargs
    )
    return xr.decode_cf(ds, use_cftime=use_cftime)


def horizontal_dims(da):
    for dim in da.dims:
        if "lon" in dim:
            lon_dim = dim
        if "lat" in dim:
            lat_dim = dim
    return (lon_dim, lat_dim)


def intersect(lamgm, phigm, lamem, phiem):
    gcm_dims = list(horizontal_dims(lamgm))
    rcm_dims = list(horizontal_dims(lamem))
    rcm_dims.append("pos")
    out_dims = rcm_dims
    # print(gcm_dims)
    # print(rcm_dims)
    # print(out_dims)
    result = xr.apply_ufunc(
        intf.intersection_points,  # first the function
        lamgm * 1.0 / 57.296,  # now arguments in the order expected by 'druint'
        phigm * 1.0 / 57.296,
        lamem * 1.0 / 57.296,
        phiem * 1.0 / 57.296,
        input_core_dims=[
            gcm_dims,
            gcm_dims,
            rcm_dims,
            rcm_dims,
        ],  # list with one entry per arg
        output_core_dims=[out_dims, out_dims],  # returned data has 3 dimensions
        #   vectorize=True,  # loop over non-core dims, in this case: time
        #   exclude_dims=set(("lev",)),  # dimensions allowed to change size. Must be a set!
        dask="parallelized",
        output_dtypes=[lamgm.dtype],
    )
    return result


def interpolate_horizontal(da, lon, lat, name=None, igr=0):
    if name is None:
        name = da.name
    lamem, phiem = lon, lat
    # lamem = xr.concat([lamem, lamem, lamem], dim='pos')
    # phiem = xr.concat([phiem, phiem, phiem], dim='pos')
    lat2d, lon2d = xr.broadcast(ds.lat, ds.lon)
    lamgm, phigm = lon2d, lat2d
    indii, indjj = intersect(lamgm, phigm, lamem, phiem)
    return interp_horiz(
        da,
        lamgm,
        phigm,
        lamem.isel(pos=igr),
        phiem.isel(pos=igr),
        indii.isel(pos=igr),
        indjj.isel(pos=igr),
        name,
    )


def interp_horiz_2d(field, lamgm, phigm, lamem, phiem, indii, indjj, name):
    """interpolates 2d global data horizontally.

    Interpolates 2d data from the global grid to the regional grid.
    """
    from intorg import intorg

    # if indii is None or indjj is None:
    #    indii, indjj = intersection_points(lamgm, phigm, lamem, phiem)
    return intorg.hiobla(field, lamgm, phigm, lamem, phiem, indii, indjj, name)


def interp_horiz(da, lamgm, phigm, lamem, phiem, indii, indjj, name, keep_attrs=False):
    """main interface"""
    gcm_dims = list(horizontal_dims(lamgm))
    rcm_dims = list(horizontal_dims(lamem))
    input_core_dims = [
        gcm_dims,
        gcm_dims,
        gcm_dims,
        rcm_dims,
        rcm_dims,
        rcm_dims,
        rcm_dims,
        [],
    ]
    result = xr.apply_ufunc(
        interp_horiz_2d,  # first the function
        da,  # now arguments in the order expected
        lamgm * 1.0 / 57.296,
        phigm * 1.0 / 57.296,
        lamem * 1.0 / 57.296,
        phiem * 1.0 / 57.296,
        indii,
        indjj,
        name,
        input_core_dims=input_core_dims,  # list with one entry per arg
        output_core_dims=[rcm_dims],  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        #  exclude_dims=set(("lev",)),  # dimensions allowed to change size. Must be a set!
        dask="parallelized",
        dask_gufunc_kwargs={"allow_rechunk": True},
        output_dtypes=[da.dtype],
    )

    result.name = name
    # result = result.to_dataset()
    if keep_attrs:
        result.attrs = da.attrs
    # result = result.transpose(..., *spatial_dims(da)[::-1])
    return result


def geopotential(fibgm, tgm, qdgm, psgm, akgm, bkgm):
    """main interface"""
    # gcm_dims = list(spatial_dims(lamgm))
    twoD_dims = list(horizontal_dims(fibgm))
    threeD_dims = list(horizontal_dims(fibgm))
    threeD_dims.append(lev_gm)
    # lev_dims.append("lev")
    # plev_dims = list(spatial_dims(da))
    # plev_dims.append("plev")
    # nlev = a.dims[0]
    input_core_dims = [
        twoD_dims,
        threeD_dims,
        threeD_dims,
        twoD_dims,
        ["lev_2"],
        ["lev_2"],
    ]
    result = xr.apply_ufunc(
        intf.geopotential,  # first the function
        fibgm,  # now arguments in the order expected
        tgm,
        qdgm,
        psgm,
        akgm,
        bkgm,
        input_core_dims=input_core_dims,  # list with one entry per arg
        #  output_core_dims=[threeD_dims],  # returned data has 3 dimensions
        output_core_dims=[twoD_dims],  # returned data has 3 dimensions
        # vectorize=True,  # loop over non-core dims, in this case: time
        # exclude_dims=set(("lev",)),  # dimensions allowed to change size. Must be a set!
        dask="parallelized",
        #  dask_gufunc_kwargs = {'allow_rechunk':True},
        output_dtypes=[fibgm.dtype],
    )
    return result


def relative_humidity(qdgm, tgm, psgm, akgm, bkgm, qwgm=None):
    """main interface"""
    if qwgm is None:
        qwgm = xr.zeros_like(qdgm)
    twoD_dims = list(horizontal_dims(qdgm))
    threeD_dims = list(horizontal_dims(qdgm)) + ["lev_gm"]
    #  print(twoD_dims)
    # threeD_dims.append("lev")
    input_core_dims = [
        threeD_dims,
        threeD_dims,
        twoD_dims,
        [akgm.dims[0]],
        [bkgm.dims[0]],
        threeD_dims,
    ]
    result = xr.apply_ufunc(
        intf.relative_humidity,  # first the function
        qdgm,  # now arguments in the order expected
        tgm,
        psgm,
        akgm,
        bkgm,
        qwgm,
        input_core_dims=input_core_dims,  # list with one entry per arg
        #  output_core_dims=[threeD_dims],  # returned data has 3 dimensions
        output_core_dims=[threeD_dims],  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        # exclude_dims=set(("lev",)),  # dimensions allowed to change size. Must be a set!
        dask="parallelized",
        #  dask_gufunc_kwargs = {'allow_rechunk':True},
        output_dtypes=[qdgm.dtype],
    )
    return result


def geo_coords(domain):
    import numpy as np

    ll_lam = domain.rlon[0]  # * 1.0/57.296
    ll_phi = domain.rlat[0]  # * 1.0/57.296
    dlam = 0.11
    dphi = 0.11
    nlam = domain.rlon.size
    nphi = domain.rlat.size
    pollam = domain.rotated_latitude_longitude.grid_north_pole_longitude
    polphi = domain.rotated_latitude_longitude.grid_north_pole_latitude
    lamem, phiem = intf.geo_coords(
        ll_lam, ll_phi, dlam, dphi, pollam, polphi, nlam + 2, nphi + 2
    )
    lamda = xr.DataArray(
        np.rad2deg(lamem),
        dims=("rlon", "rlat", "pos"),
        coords={"rlon": fibem.rlon, "rlat": fibem.rlat},
    )
    phida = xr.DataArray(
        np.rad2deg(phiem),
        dims=("rlon", "rlat", "pos"),
        coords={"rlon": fibem.rlon, "rlat": fibem.rlat},
    )
    return lamda, phida


def get_vc(ds):
    """Reads the vertical hybrid coordinate from a dataset."""
    ak_valid = ["ap_bnds", "a_bnds"]
    bk_valid = ["b_bnds"]
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
    #    if not all([ak_bnds, bk_bnds]):
    #        print('could not identify vertical coordinate, tried: {}, {}'.format(ak_valid, bk_valid))
    #        raise Exception('incomplete input dataset')
    #        ak_bnds, bk_bnds  = (ak_bnds[:1], bk_bnds[:,1])
    nlev = ak_bnds.shape[0]
    ak = np.zeros([nlev + 1], dtype=np.float64)
    bk = np.ones([nlev + 1], dtype=np.float64)
    if ds.lev.positive == "down":
        ak[:-1] = np.flip(ak_bnds[:, 1])
        bk[:-1] = np.flip(bk_bnds[:, 1])
    else:
        ak[1:] = np.flip(ak_bnds[:, 1])
        bk[1:] = np.flip(bk_bnds[:, 1])

    return xr.DataArray(ak, dims="lev_2"), xr.DataArray(bk, dims="lev_2")


def map_sst(tos, ref_ds, resample="6H", regrid=True):
    from datetime import timedelta as td
    import xesmf as xe

    # tos_res = tos
    attrs = tos.attrs
    tos_times = (ref_ds.time.min() - td(days=1), ref_ds.time.max() + td(days=1))
    tos = tos.sel(time=slice(tos_times[0], tos_times[1]))
    # return tos_res
    tos = tos.resample(time=resample).interpolate("linear").chunk({"time": 1})
    tos = tos.sel(time=ref_ds.time)
    if regrid:
        regridder = xe.Regridder(tos, ref_ds, "nearest_s2d")
        tos = regridder(tos)
        tos.attrs.update(attrs)
    return tos


def convert_units(ds):
    """convert units for use in the preprocessor"""
    try:
        if ds.sftlf.units == "%":
            print("converting sftlf units to fractional")
            attrs = ds.sftls.attrs
            ds["sftlf"] = ds.sftlf * 0.01
            attrs["units"] = 1
            ds.sftlf.attrs = attrs
    except:
        pass
    try:
        if ds.tos.units == "degC":
            print("converting tos units to K")
            attrs = ds.tos.attrs
            ds["tos"] = ds.tos + const.absolute_zero
            attrs["units"] = "K"
            ds.sftlf.attrs = attrs
    except:
        pass
    try:
        if ds.orog.units == "m":
            print("converting orog units to fractional")
            attrs = ds.orog.attrs
            ds["orog"] = ds.orog * const.grav_const
            attrs["units"] = "1"
            ds.orog.attrs = attrs
    except:
        pass
    return ds


def gfile(datasets, ref_ds, tos=None, time_range=None):
    """Creates a virtual gfile"""
    lon, lat = horizontal_dims(ref_ds)
    if time_range is None:
        time_range = ref_ds.time
    dsets = []
    for var, f in datasets.items():
        try:
            da = open_mfdataset(f, chunks={"time": 1})[var]
            da = da.sel(time=time_range)
        except:
            da = open_mfdataset(f, chunks={})[var]
        try:
            if da.lev.positive == "down":
                da = da.reindex(lev=da.lev[::-1])
        except:
            pass
        da[lon] = ref_ds[lon]
        da[lat] = ref_ds[lat]
        dsets.append(da)
    # with dask.config.set(**{'array.slicing.split_large_chunks': True}):
    ds = xr.merge(dsets, compat="override")
    if tos is not None:
        ds["tos"] = map_sst(tos, ref_ds.sel(time=time_range))
    ds["akgm"], ds["bkgm"] = get_vc(ref_ds)
    ds = ds.rename({"lev": lev_gm})
    ds = convert_units(ds)
    ds.attrs = ref_ds.attrs
    return ds


def rotate_uv(uge, vge, uvge, vuge, lamem, phiem, pollam, polphi):
    ulamem, uphiem = lamem.isel(pos=1), phiem.isel(pos=1)
    vlamem, vphiem = lamem.isel(pos=2), phiem.isel(pos=2)
    twoD_dims = list(horizontal_dims(uge))
    input_core_dims = 4 * [twoD_dims + [lev_gm]] + 4 * [twoD_dims] + 2 * [[]]
    print(input_core_dims)
    uge_rot, vge_rot = xr.apply_ufunc(
        intf.rotate_uv,  # first the function
        uge,  # now arguments in the order expected
        vge,
        uvge,
        vuge,
        ulamem * 1.0 / 57.296,
        uphiem * 1.0 / 57.296,
        vlamem * 1.0 / 57.296,
        vphiem * 1.0 / 57.296,
        pollam,
        polphi,
        input_core_dims=input_core_dims,  # list with one entry per arg
        #  output_core_dims=[threeD_dims],  # returned data has 3 dimensions
        output_core_dims=2 * [twoD_dims + [lev_gm]],  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        # exclude_dims=set(("lev",)),  # dimensions allowed to change size. Must be a set!
        dask="parallelized",
        #  dask_gufunc_kwargs = {'allow_rechunk':True},
        output_dtypes=[uge.dtype, vge.dtype],
    )
    return uge_rot, vge_rot


def pressure_correction_em(psge, tge, arfge, fibge, fibem, akgm, bkgm, kpbl):
    twoD_dims = list(horizontal_dims(psge))
    threeD_dims = list(horizontal_dims(psge)) + [lev_gm]
    input_core_dims = (
        [twoD_dims]
        + 2 * [threeD_dims]
        + 2 * [twoD_dims]
        + [[akgm.dims[0]], [bkgm.dims[0]], []]
    )
    print(input_core_dims)
    result = xr.apply_ufunc(
        intf.pressure_correction_em,  # first the function
        psge,  # now arguments in the order expected
        tge,
        arfge,
        fibge,
        fibem,
        akgm,
        bkgm,
        kpbl,
        input_core_dims=input_core_dims,  # list with one entry per arg
        output_core_dims=[twoD_dims],  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        dask="parallelized",
        output_dtypes=[psge.dtype],
    )
    return result


def interpolate_vertical(xge, psge, ps1em, akhgm, bkhgm, akhem, bkhem, varname, kpbl):
    twoD_dims = list(horizontal_dims(psge))
    threeD_dims = list(horizontal_dims(psge)) + [lev_gm]
    input_core_dims = (
        [threeD_dims]
        + 2 * [twoD_dims]
        + [[akhgm.dims[0]], [bkhgm.dims[0]], [akhem.dims[0]], [bkhem.dims[0]], [], []]
    )
    output_core_dims = [twoD_dims + [akhem.dims[0]]]
    print(output_core_dims)
    result = xr.apply_ufunc(
        intf.interp_vert,  # first the function
        xge,  # now arguments in the order expected
        psge,
        ps1em,
        akhgm,
        bkhgm,
        akhem,
        bkhem,
        varname,
        kpbl,
        input_core_dims=input_core_dims,  # list with one entry per arg
        output_core_dims=output_core_dims,  # returned data has 3 dimensions
        # exclude_dims=set(("index",)),
        vectorize=True,  # loop over non-core dims, in this case: time
        dask="parallelized",
        output_dtypes=[xge.dtype],
    )
    return result


def pressure_correction_ge(ps1em, tem, arfem, ficge, fibem, akem, bkem):
    twoD_dims = list(horizontal_dims(psge))
    threeD_dims = list(horizontal_dims(psge)) + [lev]
    input_core_dims = (
        [twoD_dims]
        + 2 * [threeD_dims]
        + 2 * [twoD_dims]
        + [[akem.dims[0]], [bkem.dims[0]]]
    )
    print(input_core_dims)
    result = xr.apply_ufunc(
        intf.pressure_correction_ge,  # first the function
        ps1em,  # now arguments in the order expected
        tem,
        arfem,
        ficge,
        fibem,
        akem,
        bkem,
        input_core_dims=input_core_dims,  # list with one entry per arg
        output_core_dims=[twoD_dims],  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        dask="parallelized",
        output_dtypes=[ps1em.dtype],
    )
    return result


def correct_uv(uem, vem, psem, akem, bkem, ll_lam, dlam, dphi):
    ulamem, uphiem = lamem.isel(pos=1), phiem.isel(pos=1)
    vlamem, vphiem = lamem.isel(pos=2), phiem.isel(pos=2)
    twoD_dims = list(horizontal_dims(uem))
    input_core_dims = (
        2 * [twoD_dims + [lev]]
        + 1 * [twoD_dims]
        + [[akem.dims[0]], [bkem.dims[0]]]
        + 3 * [[]]
    )
    print(input_core_dims)
    uge_corr, vge_corr = xr.apply_ufunc(
        intf.correct_uv,  # first the function
        uem,  # now arguments in the order expected
        vem,
        psem,
        akem,
        bkem,
        ll_lam,
        dlam,
        dphi,
        input_core_dims=input_core_dims,  # list with one entry per arg
        #  output_core_dims=[threeD_dims],  # returned data has 3 dimensions
        output_core_dims=2 * [twoD_dims + [lev]],  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        # exclude_dims=set(("lev",)),  # dimensions allowed to change size. Must be a set!
        dask="parallelized",
        #  dask_gufunc_kwargs = {'allow_rechunk':True},
        output_dtypes=[uem.dtype, vem.dtype],
    )
    return uge_corr, vge_corr
