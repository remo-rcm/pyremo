
import xarray as xr
import cordex as cx
import json

from ..core import codes


def _get_varinfo(name):
    # fails silently
    try:
        return codes.get_dict(name)
    except:
        return None


def _get_pole(ds):
    """returns the first pole we find in the dataset"""
    pol_names = ["rotated_latitude_longitude", "rotated_pole"]
    for pol in pol_names:
        if pol in ds:
            return ds[pol]
    return None


def _get_grid_definitions(CORDEX_domain, **kwargs):
    return cx.cordex_domain(CORDEX_domain, add_vertices=True, **kwargs)


def _get_cordex_pole(CORDEX_domain):
    return cx.cordex_domain(CORDEX_domain).rotated_latitude_longitude


def _set_time_units(time, units):
    time.encoding["units"] = units
    return time


def _encode_time(time):
    """encode xarray time axis into cf values

    see https://github.com/pydata/xarray/issues/4412

    """
    return xr.conventions.encode_cf_variable(time)


def _read_cmor_table(table):
    return _read_json_file(table)


def _read_json_file(filename):
    with open(filename) as f:
        data = json.load(f)
    return data


def _get_cfvarinfo(cf_varname, table):
    data = _read_cmor_table(table)
    return data['variable_entry'][cf_varname]


def _get_time_cell_method(cf_varname, table):
    return _strip_time_cell_method(_get_cfvarinfo(cf_varname, table))


def _strip_time_cell_method(cfvarinfo):
    try:
        return cfvarinfo['cell_methods'].split('time:')[1].strip()
    except:
        return None