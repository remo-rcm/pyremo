
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


