
import os
import xarray as xr


xr.set_options(keep_attrs = True)

from .core import ( const, geo_coords, geopotential, interpolate_horizontal,
                   relative_humidity, interpolate_vertical, rotate_uv, pbl_index,
                   get_akbkem, pressure_correction_em, pressure_correction_ge, 
                   correct_uv )


vcs = ['hyai', 'hybi', 'hyam', 'hybm', 'akgm', 'bkgm', 'ak', 'bk']


def get_filename(date, expid='000000', template=None):
    if template is None:
        template = 'x{}x{}.nc'
    return template.format(expid, date.strftime(format='%Y%m%d%H'))


def to_netcdf(ads, path='', expid='000000', template=None, 
              tempfiles=None, **kwargs):
    """write dataset to netcdf
    
    by default, each timestep goes into a separate output file
    
    """
    if not os.path.isdir(path):
        os.makedirs(path)
    expand_time = [var for var, da in ads.items() if 'time' in da.dims]
    if template is None:
        template = 'a{}a{}.nc'
    dates, datasets = zip(*ads.groupby("time"))
    paths = [os.path.join(path, get_filename(date, expid, template)) 
             for date in dates]
    dsets = []
    # expand time dimension only for variables not coordinates.
    for ds in datasets:
        for var, da in ds.items():
            if var in expand_time:
                ds[var] = da.expand_dims('time')
        dsets.append(ds)
    #dsets = [dset.expand_dims('time') for dset in datasets]
    xr.save_mfdataset(dsets, paths, **kwargs)
    if tempfiles is not None:
        for f in tempfiles:
            os.remove(f)
    return paths


def to_tar(files, tar_file, mode="w"):
    import tarfile
    try:
        from tqdm.notebook import tqdm
    except:
        def tqdm(x):
            return x
    tf = tarfile.open(tar_file, mode=mode)
    for f in tqdm(files, desc='creating tarfile'):
        tf.add(f, arcname=os.path.basename(f), recursive=False)
    tf.close()
    return tar_file


def remap(gds, domain_info, vc, surflib):
    """remapping workflow
    
    This function should be similar to the ones in the
    legacy fortran preprocessor intorg.
    
    """
    
    ## curvilinear coordinaetes
    fibem = surflib.FIB * const.grav_const
    lamem, phiem = geo_coords(domain_info, fibem.rlon, fibem.rlat)
    
    ## broadcast 1d global coordinates
    lat2d, lon2d = xr.broadcast(gds.lat, gds.lon)
    lamgm, phigm = lon2d, lat2d
    
    ## horizontal interpolation
    tge = interpolate_horizontal(gds.ta, lamem, phiem, 
                                 lamgm, phigm, 'T')
    psge = interpolate_horizontal(gds.ps, lamem, phiem, 
                                  lamgm, phigm, 'PS')
    uge = interpolate_horizontal(gds.ua, lamem, phiem, 
                                 lamgm, phigm,'U', 1)
    uvge = interpolate_horizontal(gds.ua, lamem, phiem, 
                                  lamgm, phigm,'U', 2)
    vge = interpolate_horizontal(gds.va, lamem, phiem, 
                                 lamgm, phigm,'V', 2)
    vuge = interpolate_horizontal(gds.va, lamem, phiem, 
                                  lamgm, phigm,'V', 1)
    qdge = interpolate_horizontal(gds.hus, lamem, phiem, 
                                  lamgm, phigm, 'QD')
    fibge = interpolate_horizontal(gds.orog, lamem, phiem, 
                                   lamgm, phigm,'FIB')
    
    ## geopotential
    if 'time' in gds.hus.dims:
        hus = gds.hus.isel(time=0)
    else:
        hus = gds.hus
    if 'time' in gds.ta.dims:
        ta = gds.ta.isel(time=0)
    else:
        ta = gds.ta
    if 'time' in gds.ps.dims:
        ps = gds.ps.isel(time=0)
    else:
        ps = gds.ps
        
    ficgm = geopotential(gds.orog, ta, hus,
                         ps, gds.akgm, gds.bkgm).squeeze(drop=True)
    ficge = interpolate_horizontal(ficgm, lamem, phiem, 
                                   lamgm, phigm, 'FIC')

    arfgm = relative_humidity(gds.hus, gds.ta, gds.ps, 
                              gds.akgm, gds.bkgm)
    arfge = interpolate_horizontal(arfgm, lamem, phiem, 
                                   lamgm, phigm, 'AREL HUM')
    
    ## wind vector rotation
    uge_rot, vge_rot = rotate_uv(uge, vge, uvge, vuge, 
                                 lamem, phiem, 
                                 domain_info['pollon'], domain_info['pollat'])
    
    ## first pressure correction
    kpbl = pbl_index(gds.akgm, gds.bkgm)
    ps1em = pressure_correction_em(psge, tge, arfge, fibge, fibem, gds.akgm, gds.bkgm, kpbl)

    ## vertical interpolation
    akhgm = 0.5 * (gds.akgm[:-1] + gds.akgm[1:])
    bkhgm = 0.5 * (gds.bkgm[:-1] + gds.bkgm[1:])
    dakgm = gds.akgm[1:] - gds.akgm[:-1]
    dbkgm = gds.bkgm[1:] - gds.bkgm[:-1]
    
    akbkem = get_akbkem(vc)
    
    tem = interpolate_vertical(tge, psge, ps1em, akhgm, bkhgm,
                               akbkem.akh, akbkem.bkh, 'T', kpbl)
    #return tem
    arfem = interpolate_vertical(arfge, psge, ps1em, akhgm, bkhgm,
                                 akbkem.akh, akbkem.bkh, 'RF', kpbl)
    
    
    ## second pressure correction and vertical interpolation of wind
    psem = pressure_correction_ge(ps1em, tem, arfem, ficge, fibem, akbkem.ak, akbkem.bk)
    psem.name = 'PS'
    
    uem = interpolate_vertical(uge_rot, psge, psem, akhgm, bkhgm, 
                           akbkem.akh, akbkem.bkh, 'U', kpbl)
    vem = interpolate_vertical(vge_rot, psge, psem, akhgm, bkhgm, 
                           akbkem.akh, akbkem.bkh, 'V', kpbl)
    
    ## wind vector correction
    philuem = domain_info['ll_lon']
    dlamem = domain_info['dlon']
    dphiem = domain_info['dlat']
    uem_corr, vem_corr = correct_uv(uem, vem, psem, akbkem.ak, 
                                    akbkem.bk, lamem, phiem, philuem,
                                    dlamem, dphiem)
    
    return xr.merge([tem, uem_corr, vem_corr, psem, arfem])


def add_soil(gfile):
    return None

def add_sst(gfile):
    return None


# variables in a-file required
#   CALL add(BOUNDARY_TABLE, 'U'     , UR       , code=131, adims=(/IE,JE,KE, 2/), leveltype=110, kake=(/1  ,KE /), ntime=2, arakawa=ARAKAWA_RIGHT)
#   CALL add(BOUNDARY_TABLE, 'V'     , VR       , code=132, adims=(/IE,JE,KE, 2/), leveltype=110, kake=(/1  ,KE /), ntime=2, arakawa=ARAKAWA_TOP)
#   CALL add(BOUNDARY_TABLE, 'T'     , TR       , code=130, adims=(/IE,JE,KE, 2/), leveltype=110, kake=(/1  ,KE /), ntime=2)
#   CALL add(BOUNDARY_TABLE, 'QD'    , QDR      , code=133, adims=(/IE,JE,KE, 2/), leveltype=110, kake=(/1  ,KE /), ntime=2)
#   CALL add(BOUNDARY_TABLE, 'QW'    , QWR      , code=153, adims=(/IE,JE,KE, 2/), leveltype=110, kake=(/1  ,KE /), ntime=2)
#   CALL add(BOUNDARY_TABLE, 'PS'    , PSR      , code=134, adims=(/IE,JE, 2/)   , leveltype=1  , ntime=2)
#   CALL add(BOUNDARY_TABLE, 'QDBL'  , QDBLR    , code=84 , adims=(/IE,JE, 2/)   , leveltype=1  , ntime=2)
#   CALL add(BOUNDARY_TABLE, 'TSW'   , TSWECHR  , code=55 , adims=(/IE,JE, 2/)   , leveltype=1  , ntime=2)

#   CALL add(BOUNDARY_TABLE, 'TSI'   , TSIECHR  , code=56 , adims=(/IE,JE, 2/)   , leveltype=1  , ntime=2)
#   CALL add(BOUNDARY_TABLE, 'SEAICE', SEAICER  , code=210, adims=(/IE,JE, 2/)   , leveltype=1  , ntime=2)
#   !
