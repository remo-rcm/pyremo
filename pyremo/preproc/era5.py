
import os
import subprocess
import tempfile

class ECMWF():

    date_fmt = "%Y-%m-%dT%H:%M:%S"
    codemap   = {130: 'ta', 134: 'ps', 131: 'ua', 132: 'va',
                 133: 'hus',  34 : 'tos',  31 :'sic',  129 : 'orog',
                 172: 'sftlf', 139: 'tsl1', 170: 'tsl2', 183: 'tsl3',
                 238: 'tsn', 236: 'tsl4', 246: 'clw', 141: 'snw', 198: 'src', 235: 'skt' ,
                 39: 'swvl1', 40: 'swvl2', 41: 'swvl2'}


    def __init__(self,df=None, time_axis='time'):
        self.cdo = Cdo(logging=True, tempdir=os.path.join(ExpVars.scratch,'python-cdo'))
        self.time_axis = time_axis
        self.df = df
        self.df['startdate'] = pd.to_datetime(self.df['startdate'])
        self.df['enddate'] = pd.to_datetime(self.df['enddate'])
        self.tmp_ds = {}

    @property
    def file_list(self):
        return self.df['path'].values

    def threeD(self, varname):
        """soil fields are handled as 2D although they have a depth dimension.
        """
        return len(self.timestep(varname).shape) == 3 and self.timestep(varname).shape[0] > 1

    @property
    def positive_down(self):
        """ERA5 data has the correct layer numbering for REMO.
        """
        return False

    @property
    def code(self):
        """returns code from dataframe, should be unique
        """
        code = self.df['code'].unique()
        if len(code) == 1:
            return code[0]
        else:
            raise Exception('code is not unique in file list')

    @property
    def ds(self):
        return Dataset(self._ref_file())

    def get_vc(self):
        """Reads the vertical hybrid coordinate from a dataset.
        """
        return self.variables['hyai'], self.variables['hybi']

   # @cached
    def _tmp_file_by_date(self, datetime, nc):
        filename = self.get_file_by_date(datetime)
        date_str = datetime.strftime(self.date_fmt)
        logging.debug('selecting {} from {}'.format(date_str, filename))
        if nc:
            return self._nc_convert( self.cdo.seldate(date_str, input=filename ) )
        else:
            return self.cdo.seldate(date_str, input=filename )

    def _ref_file(self, filename=None, nc=True):
        """the reference dataset is used to access dataset attributes.

        This is required to give the preprocessor acces to, e.g., grid
        and calendar information.
        """
        if filename is None:
            filename = self.file_list[0]
        logging.info('getting reference dataset from {}'.format(filename))
        if nc:
            return self._nc_convert( self._reftimestep( filename ) )
        else:
            return self._reftimestep( filename )

        #logging.debug('reference dataset: {}'.format(self.ref_ds.filepath()))
   # @cached
    def _reftimestep(self, filename):
        return self.cdo.seltimestep( 1, input=filename)

  #  @cached
    def _nc_convert(self, filename):
        regular = self._to_regular(filename)
        if self.code in self.codemap:
            return self.cdo.setname(self.codemap[self.code], input=regular)
        else:
            return regular

    def _gridtype(self, filename):
        griddes = self.cdo.griddes(input=filename)
        griddes = {entry.split('=')[0].strip():entry.split('=')[1].strip() for entry in griddes if '=' in entry}
        return griddes['gridtype']

    def _get_code(self, filename):
        showcode = self.cdo.showcode(input=filename)
        print(showcode)
        return showcode

  #  @cached
    def _to_regular(self, filename):
        """converts ecmwf spectral grib data to regular gaussian netcdf.

        cdo is used to convert ecmwf grid data to netcdf depending on the gridtype:
        For 'gaussian_reduced': cdo -R
            'spectral'        : cdo sp2gpl

        This follows the recommendation from the ECMWF Era5 Documentation.
        We also invert the latitudes to stick with cmor standard.

        """
        gridtype = self._gridtype(filename)
        if gridtype == 'gaussian_reduced':
            #return self.cdo.copy(options='-R -f nc', input=filename)
            gaussian = self.cdo.setgridtype('regular', options='-f nc', input=filename)
        elif gridtype == 'spectral':
            gaussian = self.cdo.sp2gpl(options='-f nc', input=filename)
        elif gridtype == 'gaussian':
            gaussian =  self.cdo.copy(options='-f nc', input=filename)
        else:
            raise Exception('unknown grid type for conversion to regular grid: {}'.format(gridtype))
        return self.cdo.invertlat(input=gaussian)

    def get_file_by_date(self, datetime):
        """search for a file in the dataframe that should contain the datetime.
        """
        df = self.df[ (datetime >= self.df['startdate']) & (datetime <= self.df['enddate']) ]
        if len(df['path'].values) == 0:
            raise Exception('no files found for {}'.format(datetime))
        elif len(df['path'].values) > 1:
            raise Exception('date selection must be unique')
        else:
            return df['path'].values[0]

    def data_by_date(self, variable, datetime, nc=True):
        return Dataset(self._tmp_file_by_date(datetime, nc)).variables[variable][0]
    
    
def _cdo_call(options='', op='', input='', output='temp'):
    if output is None:
        output = ''
    elif output == 'temp':
        output = tempfile.NamedTemporaryFile().name
    print(output)
    call = "cdo {} {} {} {}".format(options, op, input, output)
    print(call)
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


def _to_regular(filename):
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
        gaussian = _cdo_call(op='setgridtype,regular', options='-f nc -t ecmwf', input=filename)
    elif gridtype == 'spectral':
        #gaussian = Cdo(tempdir=scratch).sp2gpl(options='-f nc', input=filename)
        #gaussian = _convert_with_cdo(filename, op='sp2gpl')
        gaussian = _cdo_call(op='sp2gpl', options='-f nc -t ecmwf', input=filename)
    elif gridtype == 'gaussian':
        #gaussian =  cdo.copy(options='-f nc', input=filename)
        gaussian = _cdo_call(op='copy', options='-f nc -t ecmwf', input=filename)
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

def _sel_date(code, date, df):
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

def _get_timestep(code, date, df):
    filename = _get_file_by_date(code, date, df)
    return _cdo_call(op='seldate,{}'.format(date), input=filename)