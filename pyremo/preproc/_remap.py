import os

import cordex as cx
import numpy as np
import xarray as xr

import pyremo as pr

xr.set_options(keep_attrs=True)

from . import physics
from .core import (
    const,
    correct_uv,
    geo_coords,
    geopotential,
    get_akbkem,
    interpolate_horizontal,
    interpolate_vertical,
    pbl_index,
    pressure_correction_em,
    pressure_correction_ge,
    relative_humidity,
    rotate_uv,
)

# variables that should have a mask with fill values
fillvars = ["TSW", "SEAICE", "TSI"]

vcs = ["hyai", "hybi", "hyam", "hybm", "akgm", "bkgm", "ak", "bk"]


def get_filename(date, expid="000000", template=None):
    if template is None:
        template = "x{}x{}.nc"
    return template.format(expid, date.strftime(format="%Y%m%d%H"))


def to_netcdf(
    ads,
    path="",
    expid="000000",
    template=None,
    tempfiles=None,
    missval=1.0e20,
    **kwargs
):
    """write dataset to netcdf

    by default, each timestep goes into a separate output file

    """
    if not os.path.isdir(path):
        os.makedirs(path)
    expand_time = [var for var, da in ads.items() if "time" in da.dims]
    if template is None:
        template = "a{}a{}.nc"
    dates, datasets = zip(*ads.groupby("time"))
    paths = [os.path.join(path, get_filename(date, expid, template)) for date in dates]
    dsets = []
    # expand time dimension only for variables not coordinates.
    for ds in datasets:
        for var, da in ds.items():
            if var in expand_time:
                ds[var] = da.expand_dims("time")
            if var in fillvars:
                ds[var].encoding["_FillValue"] = missval
            else:
                ds[var].encoding["_FillValue"] = None
        dsets.append(ds)
    # dsets = [dset.expand_dims('time') for dset in datasets]
    writer = xr.save_mfdataset(dsets, paths, **kwargs)
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
    for f in tqdm(files, desc="creating tarfile"):
        tf.add(f, arcname=os.path.basename(f), recursive=False)
    tf.close()
    return tar_file


def gm_coords(gds):
    """broadcast 1d global coordinates"""
    lat2d, lon2d = xr.broadcast(gds.lat, gds.lon)
    lamgm, phigm = lon2d, lat2d
    return lamgm, phigm


def remap(gds, domain_info, vc, surflib):
    """remapping workflow

    This function should be similar to the ones in the
    legacy fortran preprocessor intorg.

    Parameters
    ----------
    gds : xarray.Dataset
        Input model dataset containing atmospheric variables for
        downscaling, including SST. The dataset must fullfil CF conventions
        containing: `ta`, `ua`, `va`, `ps`, `tos`, `orog` and `sftlf`.

    domain_info : dict
        A dictionary containing the domain information of the target domain.

    domain_info : pandas.DataFrame
        A table with the vertical coordinate coefficients `ak` and `bk`.

    surflib : xarray.Dataset
        The surface library containing the target grid land sea mask `BLA` and
        orography `FIB`.

    Returns
    -------
    Forcing Data : xarray.core.Dataset
        Dataset containing the atmospheric forcing data interpolated to the
        target domain.

    """

    ## curvilinear coordinaetes
    # remove time dimension if there is one
    fibem = surflib.FIB.squeeze(drop=True) * const.grav_const

    lamem, phiem = geo_coords(domain_info, fibem.rlon, fibem.rlat)

    ## broadcast 1d global coordinates
    lamgm, phigm = gm_coords(gds)

    ## horizontal interpolation
    tge = interpolate_horizontal(gds.ta, lamem, phiem, lamgm, phigm, "T")
    psge = interpolate_horizontal(gds.ps, lamem, phiem, lamgm, phigm, "PS")
    uge = interpolate_horizontal(gds.ua, lamem, phiem, lamgm, phigm, "U", 1)
    uvge = interpolate_horizontal(gds.ua, lamem, phiem, lamgm, phigm, "U", 2)
    vge = interpolate_horizontal(gds.va, lamem, phiem, lamgm, phigm, "V", 2)
    vuge = interpolate_horizontal(gds.va, lamem, phiem, lamgm, phigm, "V", 1)
    qdge = interpolate_horizontal(gds.hus, lamem, phiem, lamgm, phigm, "QD")
    fibge = interpolate_horizontal(gds.orog, lamem, phiem, lamgm, phigm, "FIB")

    ## geopotential
    #     if "time" in gds.hus.dims:
    #         hus = gds.hus.isel(time=0)
    #     else:
    #         hus = gds.hus
    #     if "time" in gds.ta.dims:
    #         ta = gds.ta.isel(time=0)
    #     else:
    #         ta = gds.ta
    #     if "time" in gds.ps.dims:
    #         ps = gds.ps.isel(time=0)
    #     else:
    #         ps = gds.ps

    ficgm = geopotential(
        gds.orog, gds.ta, gds.hus, gds.ps, gds.akgm, gds.bkgm
    )  # .squeeze(drop=True)

    ficge = interpolate_horizontal(ficgm, lamem, phiem, lamgm, phigm, "FIC")

    if "clw" in gds:
        # if False:
        arfgm = relative_humidity(gds.hus, gds.ta, gds.ps, gds.akgm, gds.bkgm, gds.clw)
    else:
        arfgm = relative_humidity(gds.hus, gds.ta, gds.ps, gds.akgm, gds.bkgm)
    arfge = interpolate_horizontal(arfgm, lamem, phiem, lamgm, phigm, "AREL HUM")

    ## wind vector rotation
    uge_rot, vge_rot = rotate_uv(
        uge, vge, uvge, vuge, lamem, phiem, domain_info["pollon"], domain_info["pollat"]
    )

    ## first pressure correction
    kpbl = pbl_index(gds.akgm, gds.bkgm)
    ps1em = pressure_correction_em(
        psge, tge, arfge, fibge, fibem, gds.akgm, gds.bkgm, kpbl
    )

    ## vertical interpolation
    akhgm = 0.5 * (gds.akgm[:-1] + gds.akgm[1:])
    bkhgm = 0.5 * (gds.bkgm[:-1] + gds.bkgm[1:])
    dakgm = gds.akgm[1:] - gds.akgm[:-1]
    dbkgm = gds.bkgm[1:] - gds.bkgm[:-1]

    akbkem = get_akbkem(vc)

    tem = interpolate_vertical(
        tge, psge, ps1em, akhgm, bkhgm, akbkem.akh, akbkem.bkh, "T", kpbl
    )
    # return tem
    arfem = interpolate_vertical(
        arfge, psge, ps1em, akhgm, bkhgm, akbkem.akh, akbkem.bkh, "RF", kpbl
    )

    ## second pressure correction and vertical interpolation of wind
    psem = pressure_correction_ge(ps1em, tem, arfem, ficge, fibem, akbkem.ak, akbkem.bk)
    psem.name = "PS"

    uem = interpolate_vertical(
        uge_rot, psge, psem, akhgm, bkhgm, akbkem.akh, akbkem.bkh, "U", kpbl
    )
    vem = interpolate_vertical(
        vge_rot, psge, psem, akhgm, bkhgm, akbkem.akh, akbkem.bkh, "V", kpbl
    )

    ## wind vector correction
    philuem = domain_info["ll_lon"]
    dlamem = domain_info["dlon"]
    dphiem = domain_info["dlat"]
    uem_corr, vem_corr = correct_uv(
        uem, vem, psem, akbkem.ak, akbkem.bk, lamem, phiem, philuem, dlamem, dphiem
    )

    tsw = remap_sst(
        gds.tos,
        lamem,
        phiem,
        lamgm,
        phigm,
        blagm=np.around(gds.sftlf),
        blaem=surflib.BLA.squeeze(drop=True),
    )

    # check if gcm contains seaice, else derive from sst
    if "sic" in gds:
        seaice = remap_seaice(
            gds.sic,
            lamem,
            phiem,
            lamgm,
            phigm,
            blagm=np.around(gds.sftlf),
            blaem=surflib.BLA.squeeze(drop=True),
        )
    else:
        seaice = physics.seaice(tsw)

    water_content = physics.water_content(tem, arfem, psem, akbkem.akh, akbkem.bkh)
    tsi = physics.tsi(tsw)

    ads = xr.merge(
        [tem, uem_corr, vem_corr, psem, arfem, tsw, seaice, water_content, tsi, akbkem]
    )

    grid = get_grid(domain_info)

    ads = ads.sel(rlon=grid.rlon, rlat=grid.rlat, method="nearest")
    ads["rlon"] = grid.rlon
    ads["rlat"] = grid.rlat

    ads = xr.merge([ads, grid])

    # rename for remo to recognize
    ads = ads.rename({"ak": "hyai", "bk": "hybi", "akh": "hyam", "bkh": "hybm"})

    # set global attributes
    ads.attrs = gds.attrs

    ads.attrs["history"] = "preprocessing with pyremo = {}".format(pr.__version__)

    ads = update_attrs(ads)

    # transpose to remo convention
    return ads.transpose(..., "lev", "rlat", "rlon")


def add_soil(ads):
    return None


def remap_sst(tos, lamem, phiem, lamgm, phigm, blagm, blaem):
    return interpolate_horizontal(
        tos, lamem, phiem, lamgm, phigm, "TSW", blagm=blagm, blaem=blaem
    )


def remap_seaice(sic, lamem, phiem, lamgm, phigm, blagm, blaem):
    seaice = interpolate_horizontal(
        sic, lamem, phiem, lamgm, phigm, "SEAICE", blagm=blagm, blaem=blaem
    )
    seaice = xr.where(seaice < 0.0, 0.0, seaice)
    return seaice


def get_grid(domain_info):
    return cx.create_dataset(**domain_info)


def encoding(da, missval):
    if np.isnan(da.values).any():
        return {"_FillValue": missval}
    else:
        return {"_FillValue": None}


def update_attrs(ds):
    for var, da in ds.items():
        try:
            attrs = pr.codes.get_dict(var)
            da.attrs = {}
            da.attrs["name"] = attrs["variable"]
            da.attrs["code"] = attrs["code"]
            da.attrs["description"] = attrs["description"]
            da.attrs["units"] = attrs["units"]
            # da.attrs['layer'] = attrs['layer']
            da.attrs["grid_mapping"] = "rotated_latitude_longitude"
            da.attrs["coordinates"] = "lon lat"
        except:
            pass
    return ds


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
