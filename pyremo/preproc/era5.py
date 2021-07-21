"""Cmorizer

ERA5 cmorizer is mostly based on the official ECMWF documentation of converting GRIB to NetCDF:

https://confluence.ecmwf.int/display/OIFS/How+to+convert+GRIB+to+netCDF

"""
import os
import subprocess
import tempfile
import xarray as xr
import pandas as pd

tempdir = None#os.path.join(os.environ['SCRATCH'], '.cdo_tmp')


def init_tempdir():
    global tempdir
    try:
        tempdir = os.path.join(os.environ['SCRATCH'], '.cdo_tmp')
        if not os.path.isdir(tempdir):
            os.makedirs(tempdir)
    except:
        tempdir = None
        
init_tempdir()
        
    
varmap   = {130: 'ta', 134: 'ps', 131: 'ua', 132: 'va',
             133: 'hus',  34 : 'tos',  31 :'sic',  129 : 'orog',
             172: 'sftlf', 139: 'tsl1', 170: 'tsl2', 183: 'tsl3',
             238: 'tsn', 236: 'tsl4', 246: 'clw', 141: 'snw', 198: 'src', 
             235: 'skt' , 39: 'swvl1', 40: 'swvl2', 41: 'swvl2', 
             138: 'svo', 155: 'sd'}

codemap = {var : code for code, var in varmap.items()}

show_cdo = True



    
def _cdo_call(options='', op='', input='', output='temp'):
    if output is None:
        output = ''
    elif output == 'temp':
        output = tempfile.TemporaryDirectory(dir=tempdir).name
    if isinstance(input, list):
        input = " ".join(input)
    call = "cdo {} {} {} {}".format(options, op, input, output)
    if show_cdo: print(call)
    stdout = subprocess.Popen(call, shell=True, stdout=subprocess.PIPE).stdout.read()
    if output:
        return output
    return stdout

    
def _convert_with_cdo(f,op):
    """Convert a single file into NetCDF format.
    """
    file = os.path.basename(f)
    path = os.path.dirname(f)
    output = os.path.join(path,file+'.nc')
    _cdo_call(options='-f nc', op=op, input=f, output=output)
    return output


#def _griddes(filename):
#    cdo = Cdo(tempdir=scratch)
#    griddes = cdo.griddes(input=filename)
#    griddes = {entry.split('=')[0].strip():entry.split('=')[1].strip() for entry in griddes if '=' in entry}
#    return griddes#['gridtype']


def _griddes(filename):
    griddes = _cdo_call(options='', op='griddes', input=filename, output='').decode('utf-8').split('\n')
    griddes = {entry.split('=')[0].strip():entry.split('=')[1].strip() for entry in griddes if '=' in entry}
    return griddes#['gridtype']


def _gridtype(filename):
    return _griddes(filename)['gridtype']


def _to_regular(filename, table=''):
    """converts ecmwf spectral grib data to regular gaussian netcdf.

    cdo is used to convert ecmwf grid data to netcdf depending on the gridtype:
    For 'gaussian_reduced': cdo -R
        'spectral'        : cdo sp2gpl

    This follows the recommendation from the ECMWF Era5 Documentation.
    We also invert the latitudes to stick with cmor standard.

    """
    #from cdo import Cdo
    #cdo = Cdo(tempdir=scratch)
    gridtype = _gridtype(filename)
    if gridtype == 'gaussian_reduced':
        #return self.cdo.copy(options='-R -f nc', input=filename)
        #gaussian = cdo.setgridtype('regular', options='-f nc', input=filename)
        #gaussian = _cdo_call(op='set')
        gaussian = _cdo_call(op='setgridtype,regular', options='-f nc '+table, input=filename)
    elif gridtype == 'spectral':
        #gaussian = Cdo(tempdir=scratch).sp2gpl(options='-f nc', input=filename)
        #gaussian = _convert_with_cdo(filename, op='sp2gpl')
        gaussian = _cdo_call(op='sp2gpl', options='-f nc '+table, input=filename)
    elif gridtype == 'gaussian':
        #gaussian =  cdo.copy(options='-f nc', input=filename)
        gaussian = _cdo_call(op='copy', options='-f nc '+table, input=filename)
    else:
        raise Exception('unknown grid type for conversion to regular grid: {}'.format(gridtype))
    return _cdo_call(op='invertlat', input=gaussian)
    #return gaussian


def _split_time(filename, scratch=None):
    from cdo import Cdo
    base = os.path.basename(filename)
    if scratch is None:
        scratch = os.path.join(os.environ['SCRATCH'], '.cdo-scratch')
        if not os.path.exists(scratch):
            os.mkdir(scratch)
    #output=os.path.join(scratch, base)
    return Cdo().splitsel('1,0,5', input=filename, output=os.path.join(scratch, base))
    #_cdo_call(op='splitsel,1', input=filename, output=output)
    #return output

def _seldate(code, date, df):
    filename = _get_file_by_date(code, date, df)
    return _cdo_call(op='seldate,{}'.format(date), input=filename)


def _convert_files(files, dask=False):
    """Convert many files into regular gaussian grid
    """
    results = []
    if dask:
        from dask import delayed
    else:
        def delayed(x):
            return x
    for f in files:
        if os.path.isfile(f):
            reg = delayed(_to_regular)(f)
            results.append(reg)
    return results


def to_xarray(filename, scratch=None):
    """convert an ERA5 grib file to xarray"""
    f_split = _split_time(filename)
    f_split.sort()
    
    
def _get_row_by_date(code, date, df):
    sel = df[df.code == code]#.loc[date]
    level_type = sel.level_type.unique()
    if len(level_type) > 1:
        raise Exception('non unique selection')
    level_type = level_type[0]
    if level_type == 'model_level':
        date = "{}-{}-{}".format(date[0:4], date[5:7], date[8:10])
    elif level_type == 'surface':
        date = "{}-{}-01".format(date[0:4], date[5:7])
    sel = sel.loc[date]#.path
    return sel


def _get_file_by_date(code, date, df):
    return _get_row_by_date(code, date, df).path


def _get_timestep(code, date, df, gaussian=True):
    f = _seldate(code, date, df)
    if gaussian is True:
        return _to_regular(f)
    return f

                       
def _get_timesteps(codes, dates, df, gaussian=True, delayed=False):
    from itertools import product
    if not isinstance(dates, list):
        dates = [dates]
    if not isinstance(codes, list):
        codes = [codes]
    timesteps = {code: [] for code in codes}
    if delayed is True:
        from dask import delayed
    else:
        def delayed(x):
            return x
    for code, date in product(codes, dates):
        reg = delayed(_get_timestep)(code, date, df, gaussian)
        timesteps[code].append(reg)
    return timesteps


def _to_dataarray(code, dates, df, parallel=False, use_cftime=True, 
                  chunks={'time' : 1}, **kwargs):
    timesteps = _get_timesteps(code, dates, df, delayed=parallel)
    if parallel:
        import dask
        timesteps_ = dask.compute(timesteps)[0]
    else:
        timesteps_ = timesteps
    dsets = [xr.open_mfdataset(files, use_cftime=use_cftime, chunks=chunks, **kwargs)
             for files in timesteps_.values()]
    return xr.merge(dsets)
    
                       
def _gfile(date, df):
    from . import core
    variables = ['ps']
    temp_files = [_get_timestep(codemap[var], date, df) for var in variables]
    return temp_files


def _open_catalog(url = "/pool/data/Catalogs/mistral-era5.json"):
    import intake
    cat = intake.open_esm_datastore(url)
    level_types = ['model_level', 'surface']
    return cat.search(code=list(codemap.values()), level_type=level_types, 
                      dataType='an', frequency='hourly')


def _get_catalog_df(url = "/pool/data/Catalogs/mistral-era5.json"):
    df = _open_catalog(url).df
    df['time'] = pd.to_datetime(df.validation_date)
    df = df.set_index(['time'])
    return df


def _wind(date, df):
    """compute wind from vorticity and divergence"""
    vorticity = 'svo'
    divergence = 'sd'
    vort_tmp = _get_timestep(codemap[vorticity], date, df, gaussian=False)
    div_tmp = _get_timestep(codemap[divergence], date, df, gaussian=False)
    #ds = xr.merge([vort_ds, div_ds])
    merge = _cdo_call(op='merge', input=[vort_tmp, div_tmp])
    uv = _cdo_call(op='dv2uvl', options='-f nc', input = merge)
    uv = _cdo_call(op='invertlat', input=uv)
    return uv
    u = self.cdo.selcode(131, input=uv)
    v = self.cdo.selcode(132, input=uv)
    return self.cdo.setname('ua', input = u), self.cdo.setname('va', input = v)

def _get_code(ident):
    if type(ident) == int:
        return ident
    elif ident in codemap:
        return codemap[ident]
    else:
        message = 'Unkown code or cf variable name: ' + str(ident) + \
                      '   I know: ' + str(codemap)
        raise Exception(message)

        
def _cf_rename(ds):
    for da in ds.values():
        try:
            print(da.code)
        except:
            pass

class ERA5():
    
    
    def __init__(self, catalog_url="/pool/data/Catalogs/mistral-era5.json", 
                scratch=None):
        if scratch:
            tempdir=scratch            
        self.df = _get_catalog_df(catalog_url)
        
        
    def to_xarray(self, idents, date, parallel=False, cf=True):
        """Create an xarray dataset"""
        idents = [_get_code(ident) for ident in idents]
        return _to_dataarray(idents, date, self.df,
                            parallel=parallel)
    
    
    def get_timesteps(self, idents, dates, gaussian=True, delayed=False):
        """Returns a list of cmorized files"""
        idents = [_get_code(ident) for ident in idents]
        return _get_timesteps(idents, dates, self.df, 
                              gaussian=gaussian, 
                              delayed=delayed)
    
    def wind(self, dates):
        return _wind(dates, self.df)
    
    
    
def gfile(date, parallel=False, cmorizer=None):
    if cmorizer is None:
        cmorizer = ERA5()
    variables = ['ta', 'hus', 'ps', 'sic']
    ds = cmorizer.to_xarray(variables, date, 
                            parallel=parallel)
    