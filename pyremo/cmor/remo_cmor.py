import xarray as xr
import cordex as cx
import datetime as dt
import cftime as cfdt
import json
from dateutil import relativedelta as reld
from warnings import warn

from .derived import derivator

from .utils import _get_varinfo, _get_pole, _set_time_units, _encode_time, _get_cordex_pole, _get_time_cell_method, _get_cfvarinfo, _strip_time_cell_method


try:
    import cmor
except:
    warn("no python cmor available")

from ..core import codes

xr.set_options(keep_attrs=True)

loffsets = {"3H": dt.timedelta(hours=1, minutes=30), 
            "6H": dt.timedelta(hours=3),
            "D" : dt.timedelta(hours=12)}


time_axis_names = {"point" : "time1",
                   "mean"  : "time"}

# map mip frequencies to pandas frequencies
freq_map = {
            "1hr" : "H",
            "3hr" : "3H",
            "6hr" : "6H",
            "day" : "D"}


# Y=2000

units_convert_rules = {
    "mm": (lambda x: x * 1.0 / 86400.0, "kg m-2 s-1"),
    "kg/kg": (lambda x: x, "1"),
}

#coordinate = cx.cordex_cmor_table('coordinate')
#formula_terms = cx.cordex_cmor_table('formula_terms')
#cv = cx.cordex_cmor_table('CV')

def ensure_cftime(func):
    def wrapper(date, **kwargs):
        return func(_to_cftime(date), **kwargs)

    return wrapper


def to_cftime(date, calendar="proleptic_gregorian"):
    if type(date) == dt.date:
        date = dt.datetime.combine(date, dt.time())
    elif isinstance(date, cfdt.datetime):
        # do nothing
        return date
    return cfdt.datetime(
        date.year,
        date.month,
        date.day,
        date.hour,
        date.minute,
        date.second,
        date.microsecond,
        calendar=calendar,
    )


# def _seasons_list():
#     # dummy leap year to allow input X-02-29 (leap day)
#     seasons = [('DJF', (dt.date(Y,  1,  1),  dt.date(Y,  2, 29))),
#            ('MAM', (dt.date(Y,  3, 1),  dt.date(Y,  5, 31))),
#            ('JJA', (dt.date(Y,  6, 1),  dt.date(Y,  8, 31))),
#            ('SON', (dt.date(Y,  9, 1),  dt.date(Y, 11, 30))),
#            ('DJF', (dt.date(Y, 12, 1),  dt.date(Y, 12, 31)))]
#     return seasons


# def _get_season(date):
#     """determine the meteorological season of a date"""
#     if isinstance(date, dt.datetime):
#         date = date.date()
#     date = date.replace(year=Y)
#     return next(season for season, (start, end) in _seasons_list()
#                 if start <= date <= end)


def _get_loffset(time):
    return loffsets.get(time, None)


def _clear_time_axis(ds):
    """Delete timesteps with NaN arrays"""
    for data_var in ds.data_vars:
        ds = ds.dropna(dim="time", how="all")
    return ds 

def _resample(
    ds, time, time_cell_method="point", label="left", time_offset=True, **kwargs
):
    """Resample a REMO variable."""
    if time_cell_method == "point":
        return ds.resample(time=time, label=label, **kwargs).nearest()#.interpolate("nearest")
    elif time_cell_method == "mean":
        if time_offset is True:
            loffset = _get_loffset(time)
        else:
            loffset = None
        return ds.resample(time=time, label=label, loffset=loffset, **kwargs).mean()
    else:
        raise Exception("unknown time_cell_method: {}".format(time_cell_method))


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


def _load_table(table):
    cmor.load_table(cx.cordex_cmor_table(table))


def _setup(table):
    cmor.setup(set_verbosity=cmor.CMOR_NORMAL, netcdf_file_action=cmor.CMOR_REPLACE)
    cmor.dataset_json(table)

def _get_time_axis_name(time_cell_method):
    """Get the name of the CMOR time coordinate"""
    return time_axis_names[time_cell_method]


def _define_axes(ds, table, time_cell_method=None):
    _load_table(table)
    if "time" in ds:
        if time_cell_method is None:
            warn('no time_cell_method given, assuming: point')
            time_cell_method = "point"
        time_values = _encode_time(ds.time).values
        time_axis_name = _get_time_axis_name(time_cell_method)
        cmorTime = cmor.axis(
            time_axis_name,
            coord_vals=time_values,
            cell_bounds=_get_bnds(time_values),
            units=ds.time.encoding["units"],
        )
    else:
        cmorTime = None
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


def _define_grid(ds, table, time_cell_method="point", grid_table="grids"):
    cmorTime, cmorLat, cmorLon = _define_axes(ds, table, time_cell_method=time_cell_method)
    _load_table(grid_table)

    cmorGrid = cmor.grid(
        [cmorLat, cmorLon], latitude=ds.lat.values, longitude=ds.lon.values
    )

    pole = _get_pole(ds)
    pole_dict = {
        "grid_north_pole_latitude": pole.grid_north_pole_latitude,
        "grid_north_pole_longitude": pole.grid_north_pole_longitude,
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
    if cmorTime is None:
        coords = [cmorGrid]
    else:
        coords = [cmorTime, cmorGrid]
    cmor_var = cmor.variable(da.name, da.units, coords)
    cmor.write(cmor_var, da.values)
    return cmor.close(cmor_var, file_name=file_name)


def _units_convert(da, table_file):
    """Convert units.
    
    Convert units according to the rules in units_convert_rules dict.
    Maybe metpy can do this also: https://unidata.github.io/MetPy/latest/tutorials/unit_tutorial.html
    
    """
    with open(cx.cordex_cmor_table(table_file)) as f:
        table = json.load(f)
    units = da.units
    cf_units = table["variable_entry"][da.name]["units"]
    if units != cf_units:
        warn("converting units {} to {}".format(units, cf_units))
        rule = units_convert_rules[units]
        da = rule[0](da)
        da.attrs["units"] = rule[1]
    return da

def _convert_cmor_to_resample_frequency(cmor_table):
    """Convert CMOR table name into resample frequency"""
    return resample_frequency[cmor_table]


def prepare_variable(
    ds,
    varname,
    CORDEX_domain=None,
    time_units="days since 1949-12-01T00:00:00",
    time_range=None,
    squeeze=True,
    allow_derive=False,
):
    """prepares a variable for cmorization."""
    is_ds = isinstance(ds, xr.Dataset)
    pole = _get_pole(ds)
    if pole is None:
        pole = _get_cordex_pole(CORDEX_domain)
    varinfo = _get_varinfo(varname)
    if varinfo is not None:
        remo_name = varinfo["variable"]
        cf_name = varinfo["cf_name"]
        if is_ds is True:
            var_ds = ds[remo_name]
        else:
            var_ds = ds
        var_ds.name = cf_name
    elif allow_derive is True:
        try:
            # assume it's a dataset with input variables for derivation.
            var_ds = derivator.derive(ds, varname)
        except:
            raise Exception("could not find or derive variable: {}".format(varname))
    else:
        raise Exception(
            "could not find {} in remo table, try allow_derive=True".format(varname)
        )
    # remove point coordinates, e.g, height2m
    if squeeze is True:
        var_ds = var_ds.squeeze(drop=True)
    if time_units is not None:
        var_ds["time"] = _set_time_units(ds.time, time_units)
    if CORDEX_domain is not None:
        var_ds = _crop_to_cordex_domain(var_ds, CORDEX_domain)
    #var_ds.attrs = ds.attrs
    return var_ds


def adjust_frequency(ds, cfvarinfo, input_freq=None):
    if input_freq is None and 'time' in ds.coords:
        input_freq = xr.infer_freq(ds.time)
    if input_freq is None:
        warn('could not determine frequency of input data, will assume it is correct.')
        return ds
    freq = freq_map[cfvarinfo['frequency']]
    if freq != input_freq:
        warn('resampling input data from {} to {}'.format(input_freq, freq))
        resample = _resample(ds, freq, time_cell_method=_strip_time_cell_method(cfvarinfo))
        return resample
    return ds
    


def cmorize_variable(
    ds, varname, cmor_table, dataset_table, allow_units_convert=False, 
    allow_resample=False, input_freq=None, CORDEX_domain=None, **kwargs
):
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
    allow_units_convert: bool
        Allow units to be converted if they do not agree with the
        units in the cmor table.
    resample: bool
        Allow to resample data. Handles both downsampling and upsampling.
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

    if CORDEX_domain is None:
        try:
            CORDEX_domain = ds.CORDEX_domain
        except:
            warn("could not identify CORDEX domain")
    ds_prep = prepare_variable(ds, varname, **kwargs)
    cfvarinfo = _get_cfvarinfo(varname, cmor_table)
    #time_cell_method = _get_time_cell_method(varname, cmor_table)
    ds_prep = adjust_frequency(ds_prep, cfvarinfo, input_freq)
    pole = _get_pole(ds_prep)
    if pole is None:
        pole = _get_cordex_pole(CORDEX_domain)
        ds_prep = xr.merge([ds_prep, pole])
    #return ds_prep
    if allow_units_convert is True:
        ds_prep[varname] = _units_convert(ds_prep[varname], cmor_table)
    _setup(dataset_table)
    time_cell_method = _strip_time_cell_method(cfvarinfo)
    cmorTime, cmorGrid = _define_grid(ds_prep, cmor_table, time_cell_method)
    return _cmor_write(ds_prep[varname], cmor_table, cmorTime, cmorGrid)
