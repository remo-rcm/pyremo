
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

def daily_sum(da):
    """
    Function to compute daily sums with a simple groupby approach
    """
    
    return da.groupby('time.day').sum(dim='time')

def seasonal_mean(da):
    """
    Function to compute seasonal means with a simple groupby approach
    """
    
    return da.groupby('time.season').mean(dim='time')


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
    
    
def get_regridder(finer, coarser, method='bilinear', **kwargs):
    """
    Function to regrid data bilinearly to a coarser grid
    """
    
    import xesmf as xe
    
    return xe.Regridder(finer, coarser, method=method, **kwargs)


def compare_seasons(ds1, ds2, orog1=None, orog2=None, do_height_correction=False,
                   regrid='ds1'):
    ds1 = ds1.copy()
    ds2 = ds2.copy()
    ds1_seasmean = seasonal_mean(ds1)
    ds2_seasmean = seasonal_mean(ds2)
    if regrid == "ds1":
        regridder = get_regridder(ds1, ds2)
        print(regridder)
        ds1_seasmean = regridder(ds1_seasmean)
    elif regrid == "ds2":
        regridder = get_regridder(ds2, ds1)
        print(regridder)
        ds2_seasmean = regridder(ds2_seasmean)

    if do_height_correction is True:
        orog1 = regridder(orog1)
        ds1_seasmean += height_correction(orog1, orog2)
    return ds1_seasmean - ds2_seasmean    
    #return xr.where(ds1_seasmean.mask, ds2_seasmean - ds1_seasmean, np.nan)

