import xarray as xr
import warnings

from . import _defaults as dftl

try:
    from pydruint import _druint_verip
except ModuleNotFoundError:
    warnings.warn(
        "The pressure interpolation requires installation of https://git.gerics.de/REMO/pydruint"
    )


# def prsint(ds, plevs=dftls.plevs, variables=dftl.variables,
#        t='T', ps='PS', fib='FIB', ak='hyai', bk='hybi', time='time'):
#    tda = ds[t]
#    psda = ds[ps]
#    fibda = ds[fib]
#    akda = ds[ak]
#    bkda = ds[bk]


def druint(data, t, ps, fib, plev, ak, bk, varname):
    return _druint_verip(data, t, ps, fib, ak, bk, varname, plev, lfilter=True)


def plev_coord(plev):
    plev_attrs = {
        "units": "Pa",
        "axis": "Z",
        "positive": "down",
        "long_name": "pressure",
        "standard_name": "air_pressure",
    }

    plev_coord = xr.DataArray(
        [100.0 * p for p in plev],
        dims="plev",
    )
    plev_coord.attrs = plev_attrs
    plev_coord.name = "plev"
    return plev_coord


def spatial_dims(da):
    for dim in da.dims:
        if "lon" in dim:
            lon_dim = dim
        if "lat" in dim:
            lat_dim = dim
    return (lon_dim, lat_dim)


def pressure_interpolation(da, plev, t, ps, orog, a, b, keep_attrs=False):
    """Pressure interpolation of data on 3D model levels.

    Parameters
    ----------
    da : xarray.DataArray
        Varialbe data on model levels.
    plev : array like or list
        Pressure levels to interpolate to [hPa].
    t : xarray.DataArray
        Atmospheric temperature on model levels.
    ps : xarray.DataArray
        Surface pressure.
    a : xarray.DataArray
        Hybrid sigma A coefficient at full levels or
        level interfaces.
    b : xarray.DataArrays
        Hybrid sigma B coefficient at full levels or
        level interfaces.

    Returns
    -------
    data on pressure levels : xarray.DataArray
        Returns data array interpolated to pressure levels.

    """
    srf_dims = list(spatial_dims(da))
    lev_dims = list(spatial_dims(da))
    lev_dims.append("lev")
    plev_dims = list(spatial_dims(da))
    plev_dims.append("plev")
    nlev = a.dims[0]
    result = xr.apply_ufunc(
        druint,  # first the function
        da,  # now arguments in the order expected by 'druint'
        t,
        ps,
        orog,
        plev,
        a,
        b,
        da.name,
        input_core_dims=[
            lev_dims,
            lev_dims,
            srf_dims,
            srf_dims,
            ["plev"],
            [nlev],
            [nlev],
            [],
        ],  # list with one entry per arg
        output_core_dims=[plev_dims],  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        exclude_dims=set(("lev",)),  # dimensions allowed to change size. Must be a set!
        dask="parallelized",
        output_dtypes=[da.dtype],
    )

    result.name = da.name
    # result = result.to_dataset()
    if keep_attrs:
        result.attrs = da.attrs
    result["plev"] = plev_coord(plev)
    result = result.transpose(..., "plev", *spatial_dims(da)[::-1])
    return result
