import xarray as xr

from . import core


def _horizontal_dims(da):
    for dim in da.dims:
        if "lon" in dim:
            lon_dim = dim
        if "lat" in dim:
            lat_dim = dim
    return (lon_dim, lat_dim)


def _vertical_dim(da):
    for dim in da.dims:
        if "lev" in dim:
            lev_dim = dim
        if "nhy" in dim:
            lev_dim = dim
    return lev_dim


def pressure(ps, ak, bk):
    """computes pressure at model levels.

    Uses surface pressure and vertical hybrid coordinates.

    Parameters
    ----------
    ps : xarray.DataArray
        Surface pressure field.
    ak, bk : xarray.DataArrays
        Arrays of hybrid coordinates at full levels or
        level interfaces.

    Returns
    -------
    p : xarray.DataArray
        Returns atmopspheric pressure.

    """
    return ak + bk * ps


def relative_humidity(t, qd, p, qw=None):
    """computes relative humidity from temperature, pressure and specific humidty (qd)

    algorithm from RemapToRemo addgm

    Parameters
    ----------
    t : xarray.DataArray
        temperature field [K]
    qd : xarray.DataArray
        specific humidity [kg/kg]
    p : xarray.DataArray
        atmospheric pressure [Pa]
    qw : xarray.DataArray
        liquid water content [kg/kg]

    Returns
    -------
    relhum : xarray.DataArray
        relative humidity [%]
    """
    if qw is None:
        qw = xr.zeros_like(t)
    t_dims = list(_horizontal_dims(t)) + [_vertical_dim(t)]
    p_dims = list(_horizontal_dims(p)) + [_vertical_dim(p)]
    input_core_dims = 2 * [t_dims] + [p_dims] + [t_dims]
    output_core_dims = [t_dims]
    result = xr.apply_ufunc(
        core.compute_arfgm,  # first the function
        t,  # now arguments in the order expected
        qd,
        p,
        qw,
        input_core_dims=input_core_dims,  # list with one entry per arg
        output_core_dims=output_core_dims,  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        dask="parallelized",
        output_dtypes=[t.dtype],
    )
    return result


def specific_humidity(t, relhum, p, set_meta=True):
    """computes specific humidity from temperature, pressure and relative humidity.

    algorithm from RemapToRemo addgm

    Parameters
    ----------
    t : xarray.DataArray
        temperature field [K]
    relhum : xarray.DataArray
        relative humidity [%]
    p : xarray.DataArray
        atmospheric pressure [Pa]

    Returns
    -------
    qd : xarray.DataArray
        specific humidity [%]
    """
    t_dims = list(_horizontal_dims(t)) + [_vertical_dim(t)]
    p_dims = list(_horizontal_dims(p)) + [_vertical_dim(p)]
    input_core_dims = 2 * [t_dims] + [p_dims]
    output_core_dims = [t_dims]
    print(input_core_dims)
    result = xr.apply_ufunc(
        core.specific_humidity,  # first the function
        t,  # now arguments in the order expected
        relhum,
        p,
        input_core_dims=input_core_dims,  # list with one entry per arg
        output_core_dims=output_core_dims,  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        dask="parallelized",
        output_dtypes=[t.dtype],
    )
    if set_meta is True:
        result.name = "QD"

    return result


def liquid_water_content(t, relhum, p, set_meta=True):
    """computes liquid water content from temperature, pressure and relative humidity.

    algorithm from RemapToRemo addgm

    Parameters
    ----------
    t : xarray.DataArray
        temperature field [K]
    relhum : xarray.DataArray
        relative humidity [%]
    p : xarray.DataArray
        atmospheric pressure [Pa]

    Returns
    -------
    qw : xarray.DataArray
        liquid water content [kg/kg]
    """
    t_dims = list(_horizontal_dims(t)) + [_vertical_dim(t)]
    p_dims = list(_horizontal_dims(p)) + [_vertical_dim(p)]
    input_core_dims = 2 * [t_dims] + [p_dims]
    output_core_dims = [t_dims]
    print(input_core_dims)
    result = xr.apply_ufunc(
        core.liquid_water_content,  # first the function
        t,  # now arguments in the order expected
        relhum,
        p,
        input_core_dims=input_core_dims,  # list with one entry per arg
        output_core_dims=output_core_dims,  # returned data has 3 dimensions
        vectorize=True,  # loop over non-core dims, in this case: time
        dask="parallelized",
        output_dtypes=[t.dtype],
    )
    if set_meta is True:
        result.name = "QW"

    return result
