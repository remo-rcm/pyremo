

import numpy as np
import xarray as xr

soil_temps = ['TS', 'TSL', 'TEMP2', 'TSN', 'TD3', 'TD4', 'TD5']


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


def weighted_field_mean(ds, lon='rlon', lat='rlat'):
    """
    function to compute area-weighted spatial means
    """
    weight = np.cos(np.deg2rad(ds[lat]))
    return ds.weighted(weight).mean(dim=(lon, lat))


def seasonal_mean(da):
    """
    Function to compute seasonal means with a simple groupby approach
    """
    
    return da.groupby('time.season').mean(dim='time')



class Experiment():
    
    def __init__(path):
        return None