
import os
import glob

import numpy as np
import xarray as xr

from ..archive import archive
from ..core.remo_ds import preprocess as remo_preprocess
from ..core import codes

soil_temps = ['TS', 'TSL', 'TEMP2', 'TSN', 'TD3', 'TD4', 'TD5']


def open_mfdataset(files, use_cftime=True, parallel=True, data_vars='minimal', chunks={'time':1}, 
                   coords='minimal', compat='override', drop=None, **kwargs):
    """optimized function for opening large cf datasets.

    based on https://github.com/pydata/xarray/issues/1385#issuecomment-561920115
    
    """
    def drop_all_coords(ds):
        return ds.reset_coords(drop=True)
    ds = xr.open_mfdataset(files, parallel=parallel, decode_times=False, combine='by_coords', 
                       preprocess=drop_all_coords, decode_cf=False, chunks=chunks,
                      data_vars=data_vars, coords=coords, compat=compat, **kwargs)
    return xr.decode_cf(ds, use_cftime=use_cftime)



class VariableSet():
    
    def __init__(self, varnames, name=None, long_name=None):
        self.varnames = varnames
        
    def extract(self, ds):
        res = xr.merge([ds[var] for var in self.varnames if var in ds])
        res.attrs = ds.attrs
        return res

    def stack(self, ds):
        return stack_variables(ds, self.varnames)

    
soil = VariableSet(soil_temps, name='soil_temperature', long_name='soil temperature')
        



def stack_variables(ds, varnames, name=None, long_name=None):
    found = [ds[var] for var in varnames if var in ds]
    dim = xr.DataArray(data=[var.name for var in found], dims='var', name='var')
    stack = xr.concat(found, dim=dim)
    if name is not None:
        stack.name = name
    if long_name is not None:
        stack.attrs['long_name'] = long_name
    return stack


def weighted_field_mean(ds, lon='rlon', lat='rlat', weights=None):
    """
    function to compute area-weighted spatial means
    """
    if weights is None:
        weights = np.cos(np.deg2rad(ds[lat]))
    return ds.weighted(weights).mean(dim=(lon, lat))


def seasonal_mean(da):
    """
    Function to compute seasonal means with a simple groupby approach
    """
    
    return da.groupby('time.season').mean(dim='time')


def height_correction(height1, height2):
    """returns height correction in m"""
    return (height1 - height2) * 0.0065


class Dataset():
    
    def __init__(self, mask='tas', **kwargs):
        self.filenames = {}
        self._mask_var = mask
        for key, value in kwargs.items():
            if not isinstance(value, list):
                value = [value]
            self.filenames[key] = value
        self.ds = self._init_dataset(self.filenames)
        
        
    def _init_dataset(self, filenames):
        da = []
        for var, filename in self.filenames.items():
            da.append(xr.open_mfdataset(filename)[var])
        ds = xr.merge(da).rename(self.inv_map)
        ds["mask"] = xr.where(~np.isnan(ds[self._mask_var].isel(time=0)),1,0)
        return ds
        
    def _get_id(self, id):
        return self.varmap[id]
    
    @property
    def mask(self):
        return self.ds.mask

    @property
    def inv_map(self):
        return {v: k for k, v in self.varmap.items()}
    
    @property
    def orog(self):
        return xr.open_dataset(self.filenames['orog'][0])
    
    def dataset(self, variables=None, mask=False):
        if variables is None:
            ds = self.ds
        else:
            if not isinstance(variables, list):
                variables = [variables]
            #names = [self._get_id(id) for id in variables]
            ds = xr.merge([self.ds[v] for v in variables])
        return ds
    

class CRU_TS4(Dataset):
    
    varmap = {'tas': 'tmp',
              'pr': 'pre',
              'orog': 'topo'}
    
    def __init__(self):
        path = "/mnt/lustre02/work/ch0636/eddy/pool/obs/cru/CRU/TS4.04/original"
        template = "cru_ts4.04.1901.2019.{variable}.dat.nc"
        variables = ["tmp", "pre", "cld", "dtr", "frs", "pet"]
        kwargs = {key: os.path.join(path, template.format(variable=key)) for key in variables}
        Dataset.__init__(self, 
                         **kwargs,
                         topo="/mnt/lustre02/work/ch0636/eddy/pool/obs/cru/CRU/TS4.04/original/cru404_c129.nc")

def inv_map(dict):
    return {v: k for k, v in dict.items()}
    
    
def create_dataset(filenames, mask=None, varmap=None):
    ds = xr.merge([xr.open_dataset(f)[v] for v, f in filenames.items()])
    if mask is not None:
        ds["mask"] = xr.where(~np.isnan(ds[mask].isel(time=0)),1,0)
    if varmap is not None:
        ds = ds.rename(inv_map(varmap))
    return ds
    
    
def cru_ts4():
    varmap = {'tas': 'tmp',
              'pr': 'pre',
              'orog': 'topo'}
    variables = ["tmp", "pre", "cld", "dtr", "frs", "pet"]
    path = "/mnt/lustre02/work/ch0636/eddy/pool/obs/cru/CRU/TS4.04/original"
    template = "cru_ts4.04.1901.2019.{variable}.dat.nc"
    filenames = {key: os.path.join(path, template.format(variable=key)) for key in variables}
    filenames['topo'] = "/mnt/lustre02/work/ch0636/eddy/pool/obs/cru/CRU/TS4.04/original/cru404_c129.nc"
    return create_dataset(filenames, mask='tmp', varmap=varmap)


def eobs(version="v22.0e"):
    varmap = {'tas': 'tg',
              'pr': 'rr',
              'tasmax': 'tx',
              'tasmin': 'tn',
              'rsds': 'qq',
              'psl': 'pp',
              'orog': 'elevation'}
    variables = ["tg", "tx", "tn", "rr", "qq", "pp"]
    path = "/mnt/lustre02/work/ch0636/eddy/pool/obs/eobs/{version}/original_025/day/var/{cf_name}/"
    template = "{variable}_ens_mean_0.25deg_reg_{version}.nc"
    filenames = {key: os.path.join(path, template).format(variable=key, version=version, cf_name=inv_map(varmap)[key]) 
                 for key in variables}
    filenames['elevation'] = "/mnt/lustre02/work/ch0636/eddy/pool/obs/eobs/{version}/original_025/fx/orog/elev_ens_0.25deg_reg_{version}.nc".format(version=version)
    return create_dataset(filenames, mask='tg', varmap=varmap)
    
#CRU_TS4 = Dataset(tas="/mnt/lustre02/work/ch0636/eddy/pool/obs/cru/CRU/TS4.04/original/cru_ts4.04.1901.2019.tmp.dat.nc")


class RemoExperiment():
    
    def __init__(self, path, time_range=None, preprocess=True):
        self.archive = archive.RemoArchive(path, force_parse=True)
        self.first_tfile = self.archive.first_tfile()
        
        self.ds = self._init_dataset(preprocess, time_range)
        
    def _init_dataset(self, preprocess=True, time_range=None):
        monthly = open_mfdataset(self.archive.monthly(time_range=time_range))
        if preprocess is True:
            monthly = remo_preprocess(monthly)
        lsm = xr.open_dataset(self.first_tfile).BLA.squeeze(drop=True)
        orog = xr.open_dataset(self.first_tfile).FIB.squeeze(drop=True)
        ds = xr.merge([monthly, lsm, orog])
        ds["mask"] = xr.where(lsm > 0, 1, 0)
        return ds

    @property       
    def mask(self):
        return self.ds.mask
        
    @property
    def variables(self):
        return list(self.monthly.data_vars)

    def dataset(self, variables=None):
        if variables is None:
            ds = self.ds
        else:
            if not isinstance(variables, list):
                variables = [variables]
            names = [codes.get_dict(id)['variable'] for id in variables]
            ds = xr.merge([self.ds[v] for v in names])
        ds.attrs = self.monthly.attrs
        return ds
    
    def soil(self):
        return stack_variables(self.ds, soil_temps, name="soil_temperatures", long_name="soil temperatures")
    
    def spin_up(self):
        return weighted_field_mean(self.soil())

    
    #def monthly(self, time_range=None):
    #    return self.archive.monthly(time_range=time_range)


def regrid(ds_src, ds_trg, method="bilinear", **kwargs):
    import xesmf as xe
    
    regridder = xe.Regridder(ds_src, ds_trg, method=method, **kwargs)
    print(regridder)
    return regridder(ds_src)

    
def get_regridder(finer, coarser, method='bilinear', **kwargs):
    """
    Function to regrid data bilinearly to a coarser grid
    """
    
    import xesmf as xe
    
    return xe.Regridder(finer, coarser, method=method, **kwargs)   