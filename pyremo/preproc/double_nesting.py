import os
import pyintorg.remo as intorg
from pyintorg.remo import driver
import numpy as np
import xarray as xr

dynamics = ['T', 'U', 'V', 'PS', 'QD', 'QW', 'QDBL', 'TSW', 'TSI', 'SICE']
static = ['FIB', 'BLA']
soil = ['TSL', 'TSN', 'TD3', 'TD4', 'TD5', 'TD', 'TDCL', 'SN', 'WL', 'WS', 'QDBL']


def update(ds):
    ds = ds.copy()
    ds['FIB'] = ds.FIB * 1.0 / 0.10197
    return ds
    
def horizontal_dims(da):
    for dim in da.dims:
        if "lon" in dim:
            lon_dim = dim
        if "lat" in dim:
            lat_dim = dim
    return (lon_dim, lat_dim)


def get_fshape(da):
    hdims = horizontal_dims(da)
    hshape = da[hdims[0]].size * da[hdims[1]].size
    if 'lev' in da.dims:
        return (hshape, da.lev.size)
    return (hshape)


def load_grid(ds, hm, vc):
    nlon_em = ds.dims['rlon']
    nlat_em = ds.dims['rlat']
    nlev_em = ds.dims['lev']
    nlon_hm = hm['nlon']
    nlat_hm = hm['nlat']
    nlev_hm = len(vc) - 1
    driver.load_grid(nlon_em, nlat_em, nlev_em, 
                     nlon_hm, nlat_hm, nlev_hm)
    
def load_data(ds, order='F'):
    ds = ds.copy().rename({'SEAICE': 'SICE'})
    for var in dynamics:
        load_variable(ds[var], mod=intorg.mo_nem, suffix="em", order=order)
        
def load_vc(ds, vc):
    akem = ds.hyai
    bkem = ds.hybi
    akhm = vc.ak
    bkhm = vc.bk
    driver.load_vc(akem, bkem, akhm, bkhm)
    
    
def load_options(filter_oro=False, derive_ice=True):
    driver.init_options(filter_oro, derive_ice)
    
def load_soil(ds, order='F'):
    for var in soil:
        load_variable(ds[var], mod=intorg.mo_nem, suffix="em", order=order)
    
def load_static(ds, order='F'):
    for var in static:
        load_variable(ds[var], mod=intorg.mo_nem, suffix="em", order=order)
    
def load_coordinates(ds, em, hm, order):
    keys = ['pollon', 'pollat', 'dlon', 'dlat', 'll_lon', 'll_lat']
    args_em = {key+'em' : em[key] for key in keys}
    args_hm = {key+'hm' : hm[key] for key in keys}
    driver.load_coordinates(**args_em, **args_hm)
    ds = ds.copy().rename({'PHI': 'PHI', 'RLA': 'LAM'})
    load_variable(ds.PHI* 1.0/57.296, mod=intorg.mo_nem, suffix="em", order=order)
    load_variable(ds.LAM* 1.0/57.296, mod=intorg.mo_nem, suffix="em", order=order)

    
def load_surflib(ds, order):
    ds = ds.copy()
    if 'rotated_pole' in ds.data_vars:
        ds = ds.drop('rotated_pole')
    ds = ds.rename({var:var[0:3] for var in ds.data_vars}) 
    for var in ds.data_vars:
        load_variable(ds[var], intorg.mo_nhm, suffix="hm", order=order)

        
def load_variable(da, mod, suffix, order):
    fvar = "{}{}".format(da.name.lower(), suffix)
    #setattr(mod, fvar, da.to_numpy().T.reshape(get_fshape(da), order=order).astype(np.float64))
    np.copyto(getattr(mod, fvar),  da.to_numpy().T.reshape(get_fshape(da), order=order).astype(np.float64))
    
    
def retrieve_variable(mod, name, shape, order='C'):
    array = getattr(mod, name.lower()+'hm')
    dims = ('rlat', 'rlon')
    if array.ndim == 2:
        shape = shape + (array.shape[-1],)
        dims = dims + ('lev', )
    return xr.DataArray(array.copy().reshape(shape, order=order), dims=dims, name=name)


def retrieve_forcing_data(time=None, transpose=None):
    ds = xr.merge([retrieve_variable(intorg.mo_nhm, var, (403, 363)) for var in dynamics])
    if time is not None:
        ds = ds.expand_dims(dim={'time': time}, axis=0)
    if transpose is not None:
        ds = ds.transpose(transpose)
    return ds
        
    
    
def deallocate():
    driver.deallocate_data()
    
def allocate():
    driver.allocate_data()
  
def init(ds, em, hm, vc, surflib, order='F'):
    ds = update(ds)
    surflib = update(surflib)
    load_grid(ds, hm, vc)
    allocate()
    #intorg.mo_nem.fibem[:] = ds.FIB.to_numpy().reshape(get_fshape(ds.FIB), order='F')
    load_vc(ds, vc)
    load_coordinates(ds, em, hm, order=order)
    load_static(ds, order=order)
    load_surflib(surflib, order=order)
    load_options()
    #deallocate()

    
def remap_timestep(ds):
    load_data(ds)
    load_soil(ds)
    intorg.driver.remap_remo()
    dsa = retrieve_forcing_data(time=ds.time)
    #deallocate()
    return dsa


def write_timestep(ds, path=None):
    if path is None:
        path = "./"
    filename = "a056000a{}.nc".format(ds.time.dt.strftime("%Y%m%d%H").data[0])
    filename = os.path.join(path, filename)
    ds.to_netcdf(filename)
    print('writing to: {}'.format(filename))
    return filename


def process_file(file, path=None, write=True):
    ds = pr.preprocess(xr.open_dataset(file))
    dn.init(ds, em, hm, vc, surflib)
    dsa = dn.remap_timestep(ds)
    dn.deallocate()
    if write is True:
        return dn.write_timestep(dsa, path)
    return dsa

def process_files(files, path=None, write=True):
    results = []
    for f in files:
        results.append(dask.delayed(process_file)(f, path, write))
    return results