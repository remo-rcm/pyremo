import warnings
import xarray as xr
import cordex as cx
<<<<<<< HEAD
import datetime as dt
import cftime as cfdt
=======
>>>>>>> master

try:
    import cmor
except:
    warnings.warn("no python cmor available")

from ..core import codes

xr.set_options(keep_attrs=True)

<<<<<<< HEAD
loffsets = {'3H': dt.timedelta(hours=1, minutes=30),
            '6H': dt.timedelta(hours=3)}

#Y=2000

def ensure_cftime(func):
    def wrapper(date, **kwargs):
        return func(_to_cftime(date), **kwargs)
    return wrapper


def to_cftime(date, calendar="proleptic_gregorian"):
    if type(date) == dt.date:
        print("convert")
        date = dt.datetime.combine(date, dt.time())
    elif isinstance(date, cfdt.datetime):
        # do nothing
        return date
    return cfdt.datetime(date.year, date.month, date.day, 
                         date.hour, date.minute, date.second, 
                         date.microsecond, calendar=calendar)


def _seasons_bounds(year, calendar=None):
    if calendar is None:
        import datetime as dt
        args = {}
    else:
        # calendar requires cftime
        import cftime as dt
        args = {'calendar' : calendar}
    return {'DJF': (dt.datetime(year-1, 12, 1, **args), dt.datetime(year, 3, 1, **args)),
            'MAM': (dt.datetime(year, 3, 1, **args), dt.datetime(year, 6, 1, **args)),
            'JJA': (dt.datetime(year, 6, 1, **args), dt.datetime(year, 9, 1, **args)),
            'SON': (dt.datetime(year, 9, 1, **args), dt.datetime(year, 12, 1, **args))}


def season_bounds(date):
    """Determines the temporal bounds of the meteorological season.
    
    Uses the month to determine the season and returns the
    temporal bounds of the season.
    
    Parameters
    ----------
    date : datetime object
        Date in the current season.
        
    Returns
    -------
    season : tuple
        Temporal bounds of the current meteorological season.
    
    """
    month = date.month
    if month != 12:
        year = date.year
    else:
        year = date.year + 1
    try:
        calendar = date.calendar
    except:
        calendar = None
    seasons_bounds = _seasons_bounds(year, calendar=calendar)
    return seasons_bounds[season(date)]

    
def _seasons():
    seasons = [('DJF', (12, 1, 2)),
           ('MAM', (3, 4, 5)),
           ('JJA', (6, 7, 8)),
           ('SON', (9, 10, 11))]
    return seasons


#@ensure_cftime
def season(date):
    """Determines the meteorological season.
    
    Uses the month to determine the season.
    
    Parameters
    ----------
    date : datetime object
        Date in the current season.
        
    Returns
    -------
    season : str
        Meteorological season of the current date.
    
    """
    return next(season for season, months in _seasons()
                if date.month in months)
    

def mid_of_season(date):
    """Determine the mid of the current season
    
    Parameters
    ----------
    date : datetime object
        Date in the current season.
        
    Returns
    -------
    mid_of_season : datetime object
        Mid date of the current season.
    
    """
    bounds = season_bounds(date)
    return bounds[0] + 0.5 * (bounds[1] - bounds[0])
    
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


def _resample(ds, time, time_cell_method='point', 
              label='left', time_offset=True, **kwargs):
    """Resample a REMO variable.
    """
    if time_cell_method=='point':
        return ds.resample(time=time, label=label, **kwargs).interpolate('nearest')
    elif time_cell_method=='mean':
        if time_offset is True:
            loffset = _get_loffset(time)
        else:
            loffset = None
        return ds.resample(time=time, label=label, loffset=loffset, **kwargs).mean()
    else:
        raise Exception('unknown time_cell_method: {}'.format(time_cell_method))
    
=======
>>>>>>> master

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
