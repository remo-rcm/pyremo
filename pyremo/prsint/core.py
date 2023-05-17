import warnings

import xarray as xr

try:
    from pydruint import _druint_verip
except ModuleNotFoundError:
    warnings.warn(
        "The pressure interpolation requires installation of https://gitlab.dkrz.de/remo/pydruint"
    )


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


def vertical_dim(da):
    for dim in da.dims:
        if "lev" in dim:
            return dim
    return None


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
        Variable data on model levels.
    plev : array like or list
        Pressure levels to interpolate to [hPa].
    t : xarray.DataArray
        Atmospheric temperature on model levels.
    orog : xarray.DataArray
        Orography [m].
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
    plev.sort()
    lev_dims = list(spatial_dims(da))
    lev_dims.append(vertical_dim(da))
    plev_dims = list(spatial_dims(da))
    plev_dims.append("plev")
    nlev = a.dims[0]

    t_dims = list(spatial_dims(t)) + [vertical_dim(t)]
    ps_dims = list(spatial_dims(ps))
    orog_dims = list(spatial_dims(orog))

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
            t_dims,
            ps_dims,
            orog_dims,
            ["plev"],
            [nlev],
            [nlev],
            [],
        ],  # list with one entry per arg
        output_core_dims=[plev_dims],  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        # exclude_dims=set(("lev",)),  # dimensions allowed to change size. Must be a set!
        dask="parallelized",
        output_dtypes=[da.dtype],
        keep_attrs=keep_attrs,
    )

    result.name = da.name
    result["plev"] = plev_coord(plev)
    result = result.transpose(..., "plev", *spatial_dims(da)[::-1])
    return result
