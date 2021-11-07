import warnings
import xarray as xr
import cordex as cx

try:
    import cmor
except:
    warnings.warn("no python cmor available")

from ..core import codes

xr.set_options(keep_attrs=True)


def _get_bnds(values):
    bnds = [None] * (len(values) + 1)
    bnds[0] = values[0] - (values[1] - values[0]) / 2
    bnds[len(values)] = values[-1] + (values[-1] - values[-2]) / 2
    i = 1
    while i < len(values):
        bnds[i] = values[i] - (values[i] - values[i - 1]) / 2
        i += 1
    return bnds


def _crop_to_cordex_domain(ds, domain):
    domain = cx.cordex_domain(domain)
    # the method=='nearest' approach does not work well with dask
    return ds.sel(
        rlon=slice(domain.rlon.min(), domain.rlon.max()),
        rlat=slice(domain.rlat.min(), domain.rlat.max()),
    )


def _get_varinfo(name):
    return codes.get_dict(name)


def _get_pole(ds):
    """returns the first pole we find in the dataset"""
    pol_names = ["rotated_latitude_longitude", "rotated_pole"]
    for pol in pol_names:
        return ds[pol]
    return None


def _set_time_units(time, units):
    time.encoding["units"] = units
    return time


def _encode_time(time):
    """encode xarray time axis into cf values

    see https://github.com/pydata/xarray/issues/4412

    """
    return xr.conventions.encode_cf_variable(time)


def _load_table(table):
    cmor.load_table(cx.cordex_cmor_table(table))


def _setup(table):
    cmor.setup(set_verbosity=cmor.CMOR_NORMAL, netcdf_file_action=cmor.CMOR_REPLACE)
    cmor.dataset_json(table)


def _define_axes(ds, table):
    _load_table(table)
    time_values = _encode_time(ds.time).values
    cmorTime = cmor.axis(
        "time",
        coord_vals=time_values,
        cell_bounds=_get_bnds(time_values),
        units=ds.time.encoding["units"],
    )
    cmorLat = cmor.axis(
        "gridlatitude",
        coord_vals=ds.rlat.values,
        cell_bounds=_get_bnds(ds.rlat.values),
        units=ds.rlat.units,
    )
    cmorLon = cmor.axis(
        "gridlongitude",
        coord_vals=ds.rlon.values,
        cell_bounds=_get_bnds(ds.rlon.values),
        units=ds.rlon.units,
    )
    return cmorTime, cmorLat, cmorLon


def _define_grid(ds, table, grid_table="grids"):
    cmorTime, cmorLat, cmorLon = _define_axes(ds, table)
    _load_table(grid_table)

    cmorGrid = cmor.grid(
        [cmorLat, cmorLon], latitude=ds.lat.values, longitude=ds.lon.values
    )

    pole = _get_pole(ds)
    pole_dict = {
        "grid_north_pole_latitude": 39.25,
        "grid_north_pole_longitude": -162,
        "north_pole_grid_longitude": 0.0,
    }
    cmorGM = cmor.set_grid_mapping(
        cmorGrid,
        "rotated_latitude_longitude",
        list(pole_dict.keys()),
        list(pole_dict.values()),
        ["", "", ""],
    )
    return cmorTime, cmorGrid


def _cmor_write(da, table, cmorTime, cmorGrid, file_name=True):
    cmor.load_table(cx.cordex_cmor_table(table))
    cmor_var = cmor.variable(da.name, da.units, [cmorTime, cmorGrid])
    cmor.write(cmor_var, da.values)
    return cmor.close(cmor_var, file_name=file_name)


def prepare_variable(
    ds,
    varname,
    CORDEX_domain=None,
    time_units="days since 1949-12-01T00:00:00",
    time_range=None,
    squeeze=True,
):
    """prepares a variable for cmorization."""
    if CORDEX_domain is None:
        try:
            CORDEX_domain = ds.CORDEX_domain
        except:
            warnings.warn("could not identify CORDEX domain")
    varinfo = _get_varinfo(varname)
    remo_name = varinfo["variable"]
    cf_name = varinfo["cf_name"]
    var_ds = xr.merge([ds[remo_name], _get_pole(ds)])
    var_ds = var_ds.rename_vars({remo_name: cf_name})
    # remove point coordinates, e.g, height2m
    if squeeze is True:
        var_ds = var_ds.squeeze(drop=True)
    if time_units is not None:
        var_ds["time"] = _set_time_units(ds.time, time_units)
    if CORDEX_domain is not None:
        var_ds = _crop_to_cordex_domain(var_ds, CORDEX_domain)
    var_ds.attrs = ds.attrs
    return var_ds


def cmorize_variable(ds, varname, cmor_table, dataset_table, **kwargs):
    """Cmorizes a variable.

    Parameters
    ----------
    ds : xr.Dataset
        REMO Dataset.
    varname: str
        CF name of the variable that should be cmorized.
    cmor_table : str
        Filepath to cmor table.
    dataset_table: str
        Filepath to dataset cmor table.
    **kwargs:
        Argumets passed to prepare_variable.

    Returns
    -------
    filename
        Filepath to cmorized file.


    Example
    -------
    Example for cmorization of a dataset that contains REMO output::

        $ filename = pr.cmor.cmorize_variable(ds, 'tas', 'Amon', 
                                  cx.cordex_cmor_table('remo_example'), 
                                  CORDEX_domain='EUR-11')

    """
    ds_prep = prepare_variable(ds, varname, **kwargs)
    _setup(dataset_table)
    cmorTime, cmorGrid = _define_grid(ds_prep, cmor_table)
    return _cmor_write(ds_prep[varname], cmor_table, cmorTime, cmorGrid)
