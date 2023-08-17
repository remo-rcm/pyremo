import os

import cf_xarray as cfxr  # noqa
import cordex as cx
import numpy as np
import xarray as xr

import pyremo as pr

from . import physics
from .constants import fillvars, lev_input
from .core import const, get_akbkem
from .utils import update_attrs
from .xpyintorg import (
    correct_uv,
    geo_coords,
    geopotential,
    interp_horiz_remo_cm,
    interpolate_horizontal,
    interpolate_horizontal_remo,
    interpolate_vertical,
    intersect,
    intersect_regional,
    pbl_index,
    pressure_correction_em,
    pressure_correction_ge,
    relative_humidity,
    rotate_uv,
)

# from pyremo.core.remo_ds import update_meta_info


xr.set_options(keep_attrs=True)


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
    return xr.save_mfdataset(dsets, paths, **kwargs)


def to_tar(files, tar_file, mode="w"):
    import tarfile

    try:
        from tqdm import tqdm
    except Exception:

        def tqdm(x):
            return x

    tf = tarfile.open(tar_file, mode=mode)
    for f in tqdm(files, desc="creating tarfile"):
        tf.add(f, arcname=os.path.basename(f), recursive=False)
    tf.close()
    return tar_file


def broadcast_coords(ds, coords=("lon", "lat")):
    """broadcast 1d global coordinates"""
    lat2d, lon2d = xr.broadcast(ds[coords[1]], ds[coords[0]])
    lam, phi = lon2d, lat2d
    return lam, phi


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
        Dataset containing the atmospheric and surface forcing data interpolated to the
        target domain. The dynamic fields are `T`, `U`, `V`, `QD`, `QW`, `PS`, `QDBL`,
        `TSW`, `TSI` and `SEAICE`.

    """
    # rename vertical coordinate of input to avoid conflict with output lev
    gds = gds.copy()
    gds = gds.rename({gds.cf["vertical"].name: lev_input})

    # remove time dimension if there is one
    fibem = surflib.FIB.squeeze(drop=True) * const.grav_const

    lamem, phiem = geo_coords(domain_info, fibem.rlon, fibem.rlat)

    # broadcast 1d global coordinates
    lamgm, phigm = broadcast_coords(gds)

    # compute remap matrix
    indii, indjj = intersect(lamgm, phigm, lamem, phiem)  # .compute()

    # horizontal interpolation
    tge = interpolate_horizontal(
        gds.ta, lamem, phiem, lamgm, phigm, "T", indii=indii, indjj=indjj
    )
    psge = interpolate_horizontal(
        gds.ps, lamem, phiem, lamgm, phigm, "PS", indii=indii, indjj=indjj
    )
    uge = interpolate_horizontal(
        gds.ua, lamem, phiem, lamgm, phigm, "U", 1, indii=indii, indjj=indjj
    )
    uvge = interpolate_horizontal(
        gds.ua, lamem, phiem, lamgm, phigm, "U", 2, indii=indii, indjj=indjj
    )
    vge = interpolate_horizontal(
        gds.va, lamem, phiem, lamgm, phigm, "V", 2, indii=indii, indjj=indjj
    )
    vuge = interpolate_horizontal(
        gds.va, lamem, phiem, lamgm, phigm, "V", 1, indii=indii, indjj=indjj
    )
    fibge = interpolate_horizontal(
        gds.orog, lamem, phiem, lamgm, phigm, "FIB", indii=indii, indjj=indjj
    )

    ficgm = geopotential(
        gds.orog, gds.ta, gds.hus, gds.ps, gds.akgm, gds.bkgm
    )  # .squeeze(drop=True)

    ficge = interpolate_horizontal(
        ficgm, lamem, phiem, lamgm, phigm, "FIC", indii=indii, indjj=indjj
    )

    if "clw" in gds:
        # if False:
        arfgm = relative_humidity(gds.hus, gds.ta, gds.ps, gds.akgm, gds.bkgm, gds.clw)
    else:
        arfgm = relative_humidity(gds.hus, gds.ta, gds.ps, gds.akgm, gds.bkgm)
    arfge = interpolate_horizontal(
        arfgm, lamem, phiem, lamgm, phigm, "AREL HUM", indii=indii, indjj=indjj
    )
    # return arfge
    # wind vector rotation
    uge_rot, vge_rot = rotate_uv(
        uge, vge, uvge, vuge, lamem, phiem, domain_info["pollon"], domain_info["pollat"]
    )
    # return uge_rot, vge_rot
    # first pressure correction
    kpbl = pbl_index(gds.akgm, gds.bkgm)

    ps1em = pressure_correction_em(
        psge, tge, arfge, fibge, fibem, gds.akgm, gds.bkgm, kpbl
    )

    # vertical interpolation
    akhgm = 0.5 * (gds.akgm[:-1] + gds.akgm[1:])
    bkhgm = 0.5 * (gds.bkgm[:-1] + gds.bkgm[1:])

    akbkem = get_akbkem(vc)
    ptop = akbkem.ak[0]

    tem = interpolate_vertical(
        tge, psge, ps1em, akhgm, bkhgm, akbkem.akh, akbkem.bkh, "T", kpbl, ptop=ptop
    )

    arfem = interpolate_vertical(
        arfge, psge, ps1em, akhgm, bkhgm, akbkem.akh, akbkem.bkh, "RF", kpbl, ptop=ptop
    )

    # second pressure correction and vertical interpolation of wind
    psem = pressure_correction_ge(ps1em, tem, arfem, ficge, fibem, akbkem.ak, akbkem.bk)
    psem.name = "PS"

    uem = interpolate_vertical(
        uge_rot, psge, psem, akhgm, bkhgm, akbkem.akh, akbkem.bkh, "U", kpbl, ptop=ptop
    )
    vem = interpolate_vertical(
        vge_rot, psge, psem, akhgm, bkhgm, akbkem.akh, akbkem.bkh, "V", kpbl, ptop=ptop
    )

    # correct wind with potential divergence
    philuem = domain_info["ll_lat"]
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
        blagm=xr.where(gds.tos.isnull(), 1.0, 0.0),
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
    ads.attrs["CORDEX_domain"] = domain_info.get("domain_id", "no name")

    ads = update_attrs(ads)

    # transpose to remo convention
    return ads.transpose(..., "lev", "rlat", "rlon")


def remap_remo(
    tds,
    domain_info_em,
    domain_info_hm,
    vc,
    surflib,
    initial=False,
    lice=None,
    uvcor=True,
):
    """remapping workflow for double nesting

    This function should be similar to the ones in the
    legacy fortran preprocessor intorg.

    Parameters
    ----------
    tds : xarray.Dataset
        REMO output t-file containing 3D atmospheric and 2D soil fields.
    domain_info_hm : dict
        A dictionary containing the domain information of the target domain.
    vc : pandas.DataFrame
        A table with the vertical coordinate coefficients ``ak`` and ``bk``.
    surflib : xarray.Dataset
        The surface library containing the target grid land sea mask ``BLA`` and
        orography ``FIB``.
    initial:
        If ``True``, add static and dynamic fields for initial conditions.
    uvcor: bool
        Do the u,v correction.

    Returns
    -------
    Forcing Dataset : xarray.Dataset
        Dataset containing the forcing data interpolated to the
        target domain. The dynamic variables include at least: ``T``, ``U``, ``V``, ``PS``, ``QD``,
        ``QW``, ``QDBL``, ``TSW``, ``TSI`` and ``SEAICE``.

    """
    # check if we have to derive seaice from tsw
    has_seaice = "SEAICE" in tds
    if has_seaice is True:
        tds["SEAICE"] = tds.SEAICE.fillna(0.0)
    if lice is None:
        lice = "SEAICE" not in tds
    # rename coords so they dont conflict with hm coords
    # tds = tds.copy().rename({"rlon": "rlon_em", "rlat": "rlat_em", "lev": lev_input})
    tds = tds.rename({"rlon": "rlon_em", "rlat": "rlat_em", "lev": lev_input})
    grid = get_grid(domain_info_hm)

    # curvilinear coordinaetes
    # remove time dimension if there is one
    fibhm = surflib.FIB.squeeze(drop=True) * const.grav_const
    tds["FIB"] = tds.FIB * const.grav_const

    blaem = tds.BLA
    blahm = surflib.BLA.squeeze(drop=True)
    phiem = tds.PHI.squeeze(drop=True)
    lamem = tds.RLA.squeeze(drop=True)

    lamhm, phihm = geo_coords(domain_info_hm, fibhm.rlon, fibhm.rlat)
    # lamhm = lamhm.isel(pos=0).squeeze(drop=True)
    # phihm = phihm.isel(pos=0).squeeze(drop=True)
    # return lamhm, phihm
    indemi, indemj, dxemhm, dyemhm = intersect_regional(domain_info_em, domain_info_hm)
    indemi = indemi.assign_coords(rlon=surflib.rlon, rlat=surflib.rlat)
    indemj = indemj.assign_coords(rlon=surflib.rlon, rlat=surflib.rlat)
    dxemhm = dxemhm.assign_coords(rlon=surflib.rlon, rlat=surflib.rlat)
    dyemhm = dyemhm.assign_coords(rlon=surflib.rlon, rlat=surflib.rlat)
    # return indemi, indemj, dxemhm, dyemhm

    # horizontal interpolation
    teh = interpolate_horizontal_remo(tds.T, indemi, indemj, dxemhm, dyemhm, "T")
    pseh = interpolate_horizontal_remo(tds.PS, indemi, indemj, dxemhm, dyemhm, "PS")
    pseh.name = "PSEH"
    ueh = interpolate_horizontal_remo(tds.U, indemi, indemj, dxemhm, dyemhm, "U", 1)
    # uveh = interpolate_horizontal_remo(tds.U, indemi, indemj, dxemhm, dyemhm, "U", 2)
    veh = interpolate_horizontal_remo(tds.V, indemi, indemj, dxemhm, dyemhm, "V", 4)
    # vueh = interpolate_horizontal_remo(tds.V, indemi, indemj, dxemhm, dyemhm, "V", 3)
    # qdeh = interpolate_horizontal_remo(tds.QD, indemi, indemj, dxemhm, dyemhm, "QD")
    fibeh = interpolate_horizontal_remo(tds.FIB, indemi, indemj, dxemhm, dyemhm, "FIB")

    if has_seaice is True and lice is False:
        tds["SEAICE"] = tds.SEAICE.fillna(0.0)
        siceem = tds.SEAICE
        sicehm = interp_horiz_remo_cm(
            tds.SEAICE,
            indemi,
            indemj,
            dxemhm,
            dyemhm,
            blaem,
            blahm,
            phiem,
            lamem,
            phihm,
            lamhm,
            "SICE",
        )
        sicehm.name = "SEAICE"
    else:
        # sicehm = None
        siceem = None
        sicehm = None

    tsweh = interp_horiz_remo_cm(
        tds.TSW,
        indemi,
        indemj,
        dxemhm,
        dyemhm,
        blaem,
        blahm,
        phiem,
        lamem,
        phihm,
        lamhm,
        "TSW",
        lice,
        siceem,
        sicehm,
    )

    tsieh = interp_horiz_remo_cm(
        tds.TSI,
        indemi,
        indemj,
        dxemhm,
        dyemhm,
        blaem,
        blahm,
        phiem,
        lamem,
        phihm,
        lamhm,
        "TSI",
        lice,
        siceem,
        sicehm,
    )

    if lice is True:
        sicehm = physics.seaice(tsweh)

    if initial is True:
        tds = addem_remo(tds)
        soil = remap_soil(
            tds,
            indemi,
            indemj,
            dxemhm,
            dyemhm,
            blaem,
            blahm,
            phiem,
            lamem,
            phihm,
            lamhm,
        )
    else:
        soil = xr.Dataset()

    # return pseh
    # return ueh
    ficem = geopotential(
        tds.FIB, tds.T, tds.QD, tds.PS, tds.hyai, tds.hybi
    )  # .squeeze(drop=True)

    ficeh = interpolate_horizontal_remo(ficem, indemi, indemj, dxemhm, dyemhm, "FIC")

    arfem = relative_humidity(tds.QD, tds.T, tds.PS, tds.hyai, tds.hybi, tds.QW)
    arfeh = interpolate_horizontal_remo(
        arfem, indemi, indemj, dxemhm, dyemhm, "AREL HUM"
    )

    # first pressure correction
    kpbl = pbl_index(tds.hyai, tds.hybi)

    ps1hm = pressure_correction_em(
        pseh, teh, arfeh, fibeh, fibhm, tds.hyai, tds.hybi, kpbl
    )
    ps1hm.name = "PS1HM"
    # return pseh, teh, arfeh, fibeh, fibhm
    # vertical interpolation
    akhem = tds.hyam
    bkhem = tds.hybm
    # dakem = tds.hyai[1:] - tds.hyai[:-1]
    # dbkem = tds.hybi[1:] - tds.hybi[:-1]
    akbkhm = get_akbkem(vc)
    ptop = akbkhm.ak[0]

    thm = interpolate_vertical(
        teh, pseh, ps1hm, akhem, bkhem, akbkhm.akh, akbkhm.bkh, "T", kpbl, ptop=ptop
    )
    # return thm
    # return teh,thm, pseh, ps1hm
    # return akhem, bkhem, akbkhm.akh, akbkhm.bkh

    arfhm = interpolate_vertical(
        arfeh, pseh, ps1hm, akhem, bkhem, akbkhm.akh, akbkhm.bkh, "RF", kpbl, ptop=ptop
    )

    # second pressure correction and vertical interpolation of wind
    pshm = pressure_correction_ge(ps1hm, thm, arfhm, ficeh, fibhm, akbkhm.ak, akbkhm.bk)
    pshm.name = "PS"
    # return ps1hm, pshm, thm, arfhm, ficeh, fibeh
    # return pseh, pshm
    # pseh_u = pseh.interp(rlon=pseh.rlon+0.5*0.0275, method='linear', kwargs={"fill_value": "extrapolate"})
    # pshm_u = pshm.interp(rlon=pshm.rlon+0.5*0.0275, method='linear', kwargs={"fill_value": "extrapolate"})
    # pseh_u[:,-1,:] = pseh[:,-1,:]
    # pshm_u[:,-1,:] = pshm[:,-1,:]
    uhm = interpolate_vertical(
        ueh,
        pseh.interp(
            rlon=pseh.rlon + 0.5 * domain_info_hm["dlon"],
            method="linear",
            kwargs={"fill_value": "extrapolate"},
        ),
        pshm.interp(
            rlon=pshm.rlon + 0.5 * domain_info_hm["dlon"],
            method="linear",
            kwargs={"fill_value": "extrapolate"},
        ),
        # pseh,
        # pshm,
        akhem,
        bkhem,
        akbkhm.akh,
        akbkhm.bkh,
        "U",
        kpbl,
        ptop=ptop,
    )
    # pseh_v = pseh.interp(rlat=pseh.rlat+0.5*0.0275, method='linear', kwargs={"fill_value": "extrapolate"})
    # pshm_v = pshm.interp(rlat=pshm.rlat+0.5*0.0275, method='linear', kwargs={"fill_value": "extrapolate"})
    # pseh_v[:,:,-1] = pseh[:,:,-1]
    # pshm_v[:,:,-1] = pshm[:,:,-1]
    vhm = interpolate_vertical(
        veh,
        pseh.interp(
            rlat=pseh.rlat + 0.5 * domain_info_hm["dlat"],
            method="linear",
            kwargs={"fill_value": "extrapolate"},
        ),
        pshm.interp(
            rlat=pshm.rlat + 0.5 * domain_info_hm["dlat"],
            method="linear",
            kwargs={"fill_value": "extrapolate"},
        ),
        # pseh_v,
        # pshm_v,
        akhem,
        bkhem,
        akbkhm.akh,
        akbkhm.bkh,
        "V",
        kpbl,
        ptop=ptop,
    )

    philuhm = domain_info_hm["ll_lon"]
    dlamhm = domain_info_hm["dlon"]
    dphihm = domain_info_hm["dlat"]

    if uvcor is True:
        uhm_corr, vhm_corr = correct_uv(
            uhm, vhm, pshm, akbkhm.ak, akbkhm.bk, lamhm, phihm, philuhm, dlamhm, dphihm
        )
    else:
        uhm_corr = uhm
        vhm_corr = vhm

    water_content = physics.water_content(thm, arfhm, pshm, akbkhm.akh, akbkhm.bkh)
    # tsi = physics.tsi(tsw)

    blaem = (np.around(tds.BLA),)
    blahm = surflib.BLA.squeeze(drop=True)
    phiem = tds.PHI.squeeze(drop=True)
    lamem = tds.RLA.squeeze(drop=True)

    ads = xr.merge(
        [
            thm,
            uhm_corr,
            vhm_corr,
            pshm,
            arfhm,
            water_content,
            pseh,
            tsweh,
            tsieh,
            sicehm,
            soil,
            akbkhm,
        ]
    )

    ads = ads.sel(rlon=grid.rlon, rlat=grid.rlat, method="nearest")
    ads["rlon"] = grid.rlon
    ads["rlat"] = grid.rlat

    ads = xr.merge([ads, grid])

    # rename for remo to recognize
    ads = ads.rename({"ak": "hyai", "bk": "hybi", "akh": "hyam", "bkh": "hybm"})

    if initial is True:
        ads = add_surflib(ads, surflib)
        ads = update_soil_temperatures(ads)
        ads["GLAC"] = xr.where(ads.SN > 9.1, 1.0, 0.0)

    ads = update_attrs(ads)

    # ads = update_meta_info(ads)
    # set global attributes
    ads.attrs = tds.attrs

    ads.attrs["history"] = "preprocessing with pyremo = {}".format(pr.__version__)

    # transpose to remo convention
    return ads.transpose(..., "lev", "rlat", "rlon")

    # return uhm_corr
    # return xr.merge([thm, pshm, uhm_corr, vhm_corr, qdhm, fibeh, ficeh, arfeh, ps1hm])


def remap_soil(
    tds, indemi, indemj, dxemhm, dyemhm, blaem, blahm, phiem, lamem, phihm, lamhm
):
    args = (indemi, indemj, dxemhm, dyemhm, blaem, blahm, phiem, lamem, phihm, lamhm)

    remap_vars = [
        "DTPB",
        "TSL",
        "TSN",
        "TD3",
        "TD4",
        "TD5",
        "TD",
        "TDCL",
        "WS",
        "WL",
        "SN",
    ]

    # DTPB: TEMPERATUR - DIFFERENZ TP - TB HORIZONTAL INTERPOLIEREN
    def remap_var(varname):
        return interp_horiz_remo_cm(
            tds[varname],
            *args,
            varname,
        )

    return xr.merge([remap_var(var) for var in remap_vars])


def update_soil_temperatures(ds):
    """Update land and soil temperatures in the input dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing ``T``, ``PS``, ``PSEH``, ``DTPB`` and ``TSL``

    """
    # BERECHNUNG DER SCHICHTDICKE DER PRANDTL-SCHICHT
    dpeh = ds.PSEH - pr.physics.pressure(ds.PSEH, ds.hyai[-2], ds.hybi[-2])
    dphm = ds.PS - pr.physics.pressure(ds.PS, ds.hyai[-2], ds.hybi[-2])
    ds["TSLEH"] = ds.TSL
    ds["TSL"] = ds.T.isel(lev=-1) - ds.DTPB * dphm / dpeh
    zdts = ds.TSL - ds.TSLEH
    ds["TSN"] = ds.TSN + zdts
    ds["TD3"] = ds.TD3 + zdts
    ds["TD4"] = ds.TD4 + zdts
    ds["TD5"] = ds.TD5 + zdts
    ds["TD"] = ds.TD + zdts
    ds["TDCL"] = ds.TDCL + zdts
    ds["WS"] = ds.WS * ds.WSMX
    return ds


#  dpeh(ij) = pseh(ij) - GETP(akem(KEEM),bkem(KEEM),pseh(ij),akem(1))
#  dphm(ij) = pshm(ij) - GETP(akhm(KEHM),bkhm(KEHM),pshm(ij),akhm(1))


# DO ij = 1 , IJ2HM
#   tswhm(ij) = tsweh(ij)
#   tsihm(ij) = tsieh(ij)
#   tslhm(ij) = thm(ij,KEHM) - dtpbeh(ij)*dphm(ij)/dpeh(ij)
# ENDDO
# !
# DO ij = 1 , IJ2HM
#   zdts(ij) = tslhm(ij) - tsleh(ij)
#   tsnhm(ij) = tsneh(ij) + zdts(ij)
#   td3hm(ij) = td3eh(ij) + zdts(ij)
#   td4hm(ij) = td4eh(ij) + zdts(ij)
#   td5hm(ij) = td5eh(ij) + zdts(ij)
#   tdhm(ij) = tdeh(ij) + zdts(ij)
#   tdclhm(ij) = tdcleh(ij) + zdts(ij)
# ENDDO


def addem_remo(tds):
    # ws auf relative bodenfeucht umrechnen
    tds = tds.copy()
    tds["WS"] = tds.WS.where(tds.WS < 1.0e9, 0.0) / tds.WSMX
    tds["DTPB"] = tds.T.isel({lev_input: -1}) - tds.TSL
    return tds
    # return xr.merge([wsem, dtpbem]).squeeze(drop=True)


def add_surflib(ads, surflib):
    surflib = surflib.sel(rlon=surflib.rlon[1:-1], rlat=surflib.rlat[1:-1])
    try:
        surflib = surflib.drop("rotated_pole")
    except Exception:
        pass
    return xr.merge(
        (ads, surflib.squeeze(drop=True)), join="override", compat="override"
    )


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
