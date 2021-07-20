
import os
import xarray as xr

xr.set_options(keep_attrs = True)

from .core import ( const, geo_coords, geopotential, interpolate_horizontal,
                   relative_humidity, interpolate_vertical, rotate_uv, pbl_index,
                   get_akbkem, pressure_correction_em, pressure_correction_ge, 
                   correct_uv )



def get_filename(date, expid='000000'):
    template = 'a{}a{}.nc'
    return template.format(expid, date.strftime(format='%Y%m%d%H'))


def to_netcdf(ads, path='', expid='000000'):
    """write dataset to netcdf
    
    by default, each timestep goes into a separate output file
    
    """
    dates, datasets = zip(*ads.groupby("time"))
    paths = [os.path.join(path, get_filename(date)) for date in dates]
    dsets = [dset.expand_dims('time') for dset in datasets]
    xr.save_mfdataset(dsets, paths)
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
    ficgm = geopotential(gds.orog, gds.ta.isel(time=0), gds.hus.isel(time=0),
                         gds.ps.isel(time=0), gds.akgm, gds.bkgm).squeeze(drop=True)
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
    