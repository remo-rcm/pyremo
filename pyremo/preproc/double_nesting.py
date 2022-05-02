import os

import numpy as np
import pyintorg.remo as intorg
import xarray as xr
from pyintorg.remo import driver

import pyremo as pr

from .utils import get_grid, update_attrs

dynamics = ["T", "U", "V", "PS", "QD", "QW", "QDBL", "TSW", "TSI", "SICE"]
static = ["FIB", "BLA"]
soil = ["TSL", "TSN", "TD3", "TD4", "TD5", "TD", "TDCL", "SN", "WL", "WS", "QDBL"]

output_variables = {
    "dynamics": ["T", "U", "V", "PS", "QD", "QW", "QDBL", "TSW", "TSI", "SICE"],
    "initial": [
        "TSL",
        "TSN",
        "TD3",
        "TD4",
        "TD5",
        "TD",
        "TDCL",
        "SN",
        "WL",
        "WS",
        "QDBL",
        "GLAC",
    ],
}


def update(ds):
    ds = ds.copy()
    ds["FIB"] = ds.FIB * 1.0 / 0.10197
    return ds


def horizontal_dims(da):
    for dim in da.dims:
        if "lon" in dim:
            lon_dim = dim
        if "lat" in dim:
            lat_dim = dim
    return (lon_dim, lat_dim)


def get_fshape(da):
    hdims = horizontal_dims(da)
    hshape = da[hdims[0]].size * da[hdims[1]].size
    if "lev" in da.dims:
        return (hshape, da.lev.size)
    return hshape


def load_grid(ds, hm, vc):
    nlon_em = ds.dims["rlon"]
    nlat_em = ds.dims["rlat"]
    nlev_em = ds.dims["lev"]
    nlon_hm = hm["nlon"]
    nlat_hm = hm["nlat"]
    nlev_hm = len(vc) - 1
    driver.load_grid(nlon_em, nlat_em, nlev_em, nlon_hm, nlat_hm, nlev_hm)


def load_data(ds, order="F"):
    ds = ds.copy().rename({"SEAICE": "SICE"})
    for var in dynamics:
        load_variable(ds[var], mod=intorg.mo_nem, suffix="em", order=order)


def load_vc(ds, vc):
    akem = ds.hyai
    bkem = ds.hybi
    akhm = vc.ak
    bkhm = vc.bk
    driver.load_vc(akem, bkem, akhm, bkhm)


def load_options(filter_oro=False, derive_ice=True):
    driver.init_options(filter_oro, derive_ice)


def load_soil(ds, order="F"):
    for var in soil:
        load_variable(ds[var], mod=intorg.mo_nem, suffix="em", order=order)


def load_static(ds, order="F"):
    for var in static:
        load_variable(ds[var], mod=intorg.mo_nem, suffix="em", order=order)


def load_coordinates(ds, em, hm, order):
    keys = ["pollon", "pollat", "dlon", "dlat", "ll_lon", "ll_lat"]
    args_em = {key + "em": em[key] for key in keys}
    args_hm = {key + "hm": hm[key] for key in keys}
    driver.load_coordinates(**args_em, **args_hm)
    ds = ds.copy().rename({"PHI": "PHI", "RLA": "LAM"})
    load_variable(ds.PHI * 1.0 / 57.296, mod=intorg.mo_nem, suffix="em", order=order)
    load_variable(ds.LAM * 1.0 / 57.296, mod=intorg.mo_nem, suffix="em", order=order)


def load_surflib(ds, order):
    ds = ds.copy()
    if "rotated_pole" in ds.data_vars:
        ds = ds.drop("rotated_pole")
    ds = ds.rename({var: var[0:3] for var in ds.data_vars})
    for var in ds.data_vars:
        load_variable(ds[var], intorg.mo_nhm, suffix="hm", order=order)


def load_variable(da, mod, suffix, order):
    fvar = "{}{}".format(da.name.lower(), suffix)
    np.copyto(
        getattr(mod, fvar),
        da.to_numpy().T.reshape(get_fshape(da), order=order).astype(np.float64),
    )


def retrieve_variable(mod, name, shape, order="C"):
    array = getattr(mod, name.lower() + "hm")
    dims = ("rlat", "rlon")
    if array.ndim == 2:
        shape = shape + (array.shape[-1],)
        dims = dims + ("lev",)
    return xr.DataArray(array.copy().reshape(shape, order=order), dims=dims, name=name)


def retrieve_variables(variables, dims, time=None, transpose=None, crop=True):
    ds = xr.merge([retrieve_variable(intorg.mo_nhm, var, dims) for var in variables])
    if time is not None:
        ds = ds.expand_dims(dim={"time": time}, axis=0)
    if transpose is not None:
        ds = ds.transpose(transpose)
    if crop is True:
        ds = ds.sel(rlon=ds.rlon[1:-1], rlat=ds.rlat[1:-1])
    if "SICE" in ds:
        ds = ds.rename({"SICE": "SEAICE"})
    return ds


def deallocate():
    driver.deallocate_data()


def allocate():
    driver.allocate_data()


def init(ds, em, hm, vc, surflib, order="F"):
    ds = update(ds)
    surflib = update(surflib)
    load_grid(ds, hm, vc)
    allocate()
    # intorg.mo_nem.fibem[:] = ds.FIB.to_numpy().reshape(get_fshape(ds.FIB), order='F')
    load_vc(ds, vc)
    load_coordinates(ds, em, hm, order=order)
    load_static(ds, order=order)
    load_surflib(surflib, order=order)
    load_options()


def remap_timestep(ds, hm, initial=False):
    variables = output_variables["dynamics"]
    if initial is True:
        variables += output_variables["initial"]
    load_data(ds)
    load_soil(ds)
    intorg.driver.remap_remo()
    dsa = retrieve_variables(
        variables, dims=(hm["nlat"] + 2, hm["nlon"] + 2), time=ds.time
    )
    return dsa


def write_timestep(ds, expid=None, path=None):
    if path is None:
        path = "./"
    if expid is None:
        expid = 0
    filename = "a{expid:06d}a{date}.nc".format(
        expid=int(expid), date=ds.time.dt.strftime("%Y%m%d%H").data[0]
    )
    filename = os.path.join(path, filename)
    ds.to_netcdf(filename)
    print("writing to: {}".format(filename))
    return filename


def process_file(
    file, em, hm, vc, surflib, initial=False, lice=False, path=None, write=True
):
    if isinstance(file, str):
        ds = pr.preprocess(xr.open_dataset(file))
    else:
        ds = file
    if isinstance(surflib, str):
        surflib = pr.preprocess(xr.open_dataset(surflib))

    init(ds, em, hm, vc, surflib)
    dsa = remap_timestep(ds, hm, initial)
    deallocate()
    dsa = update_dataset(dsa, ds, hm)
    if initial is True:
        # dsa = xr.merge([dsa, surflib.sel(rlon=surflib.rlon[1:-1], rlat=surflib.rlat[1:-1])])
        dsa = dsa.merge(surflib.sel(rlon=surflib.rlon[1:-1], rlat=surflib.rlat[1:-1]))
    if write is True:
        return write_timestep(dsa, path)
    return dsa


def process_files(
    files,
    em,
    hm,
    vc,
    surflib,
    initial=False,
    lice=False,
    path=None,
    write=True,
    parallel=False,
    compute=True,
):
    results = []
    if parallel is True:
        import dask

        for f in files:
            results.append(
                dask.delayed(process_file)(
                    f, em, hm, vc, surflib, initial, lice, path, write
                )
            )
        if compute is True:
            return dask.compute(*results)
        return results
    else:
        for f in files:
            results.append(
                process_file(f, em, hm, vc, surflib, initial, lice, path, write)
            )
        return results


def update_dataset(dsa, ds, hm):
    """add grid and meta data"""
    dsa = dsa.update(get_grid(hm))
    dsa = update_attrs(dsa)
    for key, value in ds.attrs.items():
        dsa.attrs[key + "_forcing"] = value
    dsa.attrs["CORDEX_domain"] = hm.get("short_name", "UNKNOWN")
    dsa.attrs["domain"] = hm.get("long_name", "UNKNOWN")
    return dsa


def remap(
    files,
    em,
    hm,
    vc,
    surflib,
    initial=False,
    lice=None,
    path=None,
    write=True,
    parallel=False,
    compute=True,
):
    """remapping workflow for double nesting

    This function should be similar to the ones in the
    legacy fortran preprocessor intorg.

    Parameters
    ----------
    files : list of filenames, xarray.Dataset
        REMO output t-file data containing 3D atmospheric and 2D soil fields.

    em : dict
        A dictionary containing the domain information of the source domain.

    hm : dict
        A dictionary containing the domain information of the target domain.

    vc : pandas.DataFrame
        A table with the vertical coordinate coefficients ``ak`` and ``bk``.

    surflib : filename or xarray.Dataset
        The surface library containing the target grid land sea mask ``BLA`` and
        orography ``FIB``.

    initial:
        If ``True``, add static and dynamic fields for initial conditions.

    Returns
    -------
    Forcing Data : filename or xarray.Dataset
        Dataset containing the forcing data interpolated to the
        target domain. The dynamic variables include at least: ``T``, ``U``, ``V``, ``PS``, ``QD``,
        ``QW``, ``QDBL``, ``TSW``, ``TSI`` and ``SEAICE``.

    """
    if not isinstance(files, (list, tuple)):
        files = [files]
    return process_files(
        files, em, hm, vc, surflib, initial, lice, path, write, parallel, compute
    )
