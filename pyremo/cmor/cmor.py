
import xarray as xr
from ..core import codes

#time_units = "days since 1949-12-01T00:00:00Z"


def _get_varinfo(name):
    return codes.get_dict(name)


def _get_pole(ds):
    """returns the first pole we find in the dataset"""
    pol_names = ['rotated_latitude_longitude', 'rotated_pole']
    for pol in pol_names:
        return ds[pol]
    return None


def _set_time_units(time, units):
    time.encoding['units'] = units
    return time


def _encode_time(time):
    """encode xarray time axis into cf values
    
    see https://github.com/pydata/xarray/issues/4412
    
    """
    return xr.conventions.encode_cf_variable(time)


def cmorize_variable(ds, varname, time_units="days since 1949-12-01T00:00:00", squeeze=True):
    """cmorize a single variable
    """
    varinfo = _get_varinfo(varname)
    remo_name = varinfo['variable']
    cf_name = varinfo['cf_name']
    var_ds = xr.merge([ds[remo_name], _get_pole(ds)])
    var_ds = var_ds.rename_vars({remo_name: cf_name})
    # remove point coordinates, e.g, height2m
    if squeeze is True:
        var_ds = var_ds.squeeze(drop=True)
    if time_units is not None:
        var_ds["time"] = _set_time_units(ds.time, time_units)
    return var_ds
