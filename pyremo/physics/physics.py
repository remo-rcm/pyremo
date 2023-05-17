import numpy as np
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


def pressure(ps, ak, bk, ptop=0.0, z_coord=None):
    """computes pressure at model levels.

    Uses surface pressure and vertical hybrid coordinates.

    Parameters
    ----------
    ps : xarray.DataArray
        Surface pressure field.
    ak : xarray.DataArray
        Hybrid sigma A coefficient at full levels or
        level interfaces.
    bk : xarray.DataArrays
        Hybrid sigma B coefficient at full levels or
        level interfaces.
    ptop : float
        Pressure at the top of the atmosphere. Defaults to
        zero.
    z_coord : xr.DataArray
        If provided, the Z coordinate will be replace by z_coord.
        Useful, if ak and bk have no coordinates but only dims.

    Returns
    -------
    p : xarray.DataArray
        Returns atmopspheric pressure.

    """
    p = ak + bk * (ps - ptop)
    if z_coord is not None:
        z_dim = (set(p.dims) - set(p.coords)).pop()
        p = p.rename({z_dim: z_coord.name}).assign_coords({z_coord.name: z_coord})
    return p


def relative_humidity(t, qd, p, qw=None, set_meta=True):
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
        # core.relative_humidity,  # first the function
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
    if set_meta:
        result.name = "RELHUM"
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


def precipitation_flux(aprl, aprc):
    """Precipitation flux ``pr`` [mm]"""
    return aprl + aprc


# remo: dew = dew2 (168)
def water_vapour(t):
    """Partial pressure of water vapour e [Pa].

    Computes water vapour from  temperature.

    references: https://doi.org/10.1175/1520-0450(1967)006%3C0203:OTCOSV%3E2.0.CO;2

    compare to: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_vapor_pressure.html#metpy.calc.saturation_vapor_pressure

    """
    T_0 = 273.15
    T_rw = 35.86  # over water
    a = 17.269
    # cdo -mulc,610.78 -exp -div -mulc,17.5 -subc,273.15 a
    return 610.78 * np.exp(a * (t - T_0) / (t - T_rw))


def specific_humidity_from_dewpoint(dew, ps):
    """Specific humidity ``huss`` [kg/kg].

    Computes specific humidity `huss` from dewpoint temperature and pressure.

    reference: https://de.wikipedia.org/wiki/Luftfeuchtigkeit#Spezifische_Luftfeuchtigkeit

    compare to: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.specific_humidity_from_dewpoint.html

    """
    e = water_vapour(dew)
    return (0.622 * e) / (ps - 0.378 * e)


def relative_humidity_from_dewpoint(dew, t2m):
    """Relative humidity ``hurs`` [%].

    Computes relative humidity ``hurs`` from dewpoint and air temperature.

    compare to: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_dewpoint.html
    """
    e_dew = water_vapour(dew)
    e_t2m = water_vapour(t2m)
    return e_dew / e_t2m


def surface_runoff_flux(runoff, drain):
    """Surface runoff ``mrros`` [mm].

    Computes surface runoff flux ``mrros`` from total runoff and drainage.
    """
    return runoff - drain


def surface_downwelling_shortwave_flux_in_air(srads, sradsu):
    """Surface downwelling shortwave flux in air ``rsds`` [W m-2].

    Computes surface downwelling shortwave flux in air ``rsds`` from net surface solar radiation and surface solar radiation upward.
    """
    return srads - sradsu


def surface_downwelling_longwave_flux_in_air(trads, tradsu):
    """Surface downwelling longwave flux in air ``rlds`` [W m-2].

    Computes surface downwelling longwave flux in air ``rlds`` from net surface thermal radiation and surface thermal radiation upward.
    """
    return trads - tradsu


def toa_incoming_shortwave_flux(srad0, srad0u):
    """TOA incoming shortwave flux ``rsdt`` [W m-2].

    Computes TOA incoming shortwave flux ``rsdt`` from net top solar radiation and top solar radiation upward.
    """
    return srad0 - srad0u


def water_evapotranspiration_flux(evap):
    """Water evapotranspiration flux ``evspsbl`` [mm].

    Computes water evapotranspiration flux ``evspsbl`` from surface evaporation.
    """
    return evap * (-1)
