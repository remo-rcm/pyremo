
import os

import pandas as pd
from pathlib import Path

import tarfile

import parse

file_pattern = "e{usr_nr:3d}{exp_nr:3d}{type:1}{date}"
efile_pattern = "e{usr_nr:3d}{exp_nr:3d}{type:1}_c{code:3d}_{date}"

date_patterns = ["%Y", "%Y%m", "%Y%m%d", "%Y%m%d", "%Y%m%d%H", "%Y%m%d%H%M"]

patterns = [efile_pattern, file_pattern]

try:
    from tqdm import tqdm
except:
    def tqdm(x):
        return x

def cdo_call(options='', op='', input='', output=None):
    call = "cdo {} {} {} {}".format(options, op, input, output)
    subprocess.Popen(call, shell=True, stdout=subprocess.PIPE).stdout.read()
    return output
    #return subprocess.run("cdo {} {} {} {}".format(options, op, input, output), shell=True)

    
def convert_with_cdo(f):
    """Convert a single file into NetCDF format.
    """
    file = os.path.basename(f)
    path = os.path.dirname(f)
    #cdo = Cdo()
    #return cdo.copy(options='-f nc', input=f, output=os.path.join(path,'nc',file+'.nc'))
    return cdo_call(options='-f nc', op='copy', input=f, output=os.path.join(path,'nc',file+'.nc'))


def convert_files(files, with_dask=False):
    """Convert many files into NetCDF format.
    """
    results = []
    if with_dask:
        from dask import delayed
    else:
        def delayed(x):
            return x
    for f in files:
        if os.path.isfile(f):
            nc = delayed(convert_with_cdo)(f)
            results.append(nc)
    return results


def get_filename_from_archive(tar, pattern):
    return next(f for f in tar.getnames() if pattern in f)


def extract_file(filename, pattern, path):
    open_tar = tarfile.open(filename, 'r')
    filename = get_filename_from_archive(open_tar, pattern)
    print(filename)
    #dest = os.path.join(scratch, filename)
    open_tar.extract(filename, path=path)
    return os.path.join(path, filename)


def extract_files(tars, pattern="", path=None, parallel=False):
    if path is None:
        path = os.path.join(os.environ['SCRATCH'], '_archive_extract')
    else:
        path = ""
    filenames = []
    if parallel is True:
        from dask import delayed
        futures = []
    else:
        def delayed(x):
            return x
        futures = None
    for tar in tars:
        filenames.append(delayed(extract_file)(tar, pattern, path))
    return filenames


def _find_files(dir, file_template="e*", exclude=None):
    """search a directory for remo output."""
    pattern = file_template
    return [path.absolute() for path in Path(dir).rglob(pattern)]


def _parse_file(file):
    for pattern in patterns:
        finfo = parse.parse(pattern, file)
        if finfo:
            return finfo.named


def get_data_catalog(dir, parse_dates=False, exclude=None):
    filepathes = _find_files(dir, exclude=exclude)
    df_all = pd.DataFrame()
    for path in tqdm(filepathes):
        finfo = _parse_file(path.stem)
        if finfo is None:
            print("could not parse: ", path.stem)
            continue
        finfo["path"] = str(path)
        finfo["suffix"] = path.suffix
        df = pd.DataFrame(finfo, index=[0])
        if parse_dates:
            for dpat in date_patterns:
                try:
                    df["date"] = pd.to_datetime(df.date, format=dpat)
                    break
                except:
                    pass
        df_all = df_all.append(df, ignore_index=True)
    try:
        df_all["code"] = df_all.code.astype(pd.Int64Dtype())
    except:
        pass
    return df_all


def _time_range(df, time_range):
    return df[(df['date'] >= time_range[0]) & (df['date'] <= time_range[1])]


def _search(df, time_range=None, **kwargs):
    for key, value in kwargs.items():
        df = df[df[key] == value]
    if time_range is not None:
            return _time_range(df, time_range)
    return df


class RemoCatalog():
    
    def __init__(self, path, force_parse=False):
        self.path = path
        self.csv = os.path.join(self.path, 'data_catalog.csv')
        if os.path.isfile(self.csv) and force_parse is False:
            self.df = pd.read_csv(self.csv)
        else:
            self.df = get_data_catalog(path, parse_dates=True)
            self.to_csv()
        
    def to_csv(self, path=None):
        self.df.to_csv(self.csv, index=False)
  
    def search(self, time_range=None, **kwargs):
        return _search(self.df, time_range=time_range, **kwargs)#.sort_values('date') 

    def archive(self, type, suffix='.tar', time_range=None, **kwargs):
        return _search(self.df, type=type, suffix=suffix,
                       time_range=time_range, **kwargs).sort_values('date')
    
    
class RemoArchive():
    
    def __init__(self, path, force_parse=False):
        self.path = path
        self.catalog = RemoCatalog(path, force_parse)
        
    def extract(self, type, code=None, time=None):
        pass
                 
    def _extract_code(self, code, time_range=None, **kwargs):
        pattern = "c{:03d}".format(code)
        tars = list(self.catalog.archive('e', time_range=time_range, **kwargs).path)
        return extract_files(tars, pattern=pattern, path=None, parallel=True)
            

    
    
        
        