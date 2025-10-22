import cf_xarray as cfxr  # noqa
import cordex as cx
import numpy as np
import xarray as xr

import pyremo as pr

from . import physics
from .constants import lev_input
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

# required variables for initial conditions:

# ['T', 'U', 'V', 'PS', 'RF', 'QW', 'QD', 'QDBL', 'PSEH', 'TSW', 'TSI', 'SEAICE', 'DTPB',
#  'TSL', 'TSN', 'TD3', 'TD4', 'TD5', 'TD', 'TDCL', 'WS', 'WL', 'SN',
#  'FIB', 'BLA', 'AZ0', 'ALB', 'VGRAT', 'VAROR', 'VLT',
#  'FOREST', 'FAO', 'WSMX', 'BETA', 'WMINLOK', 'WMAXLOK', 'TSLEH', 'GLAC']

xr.set_options(keep_attrs=True)


def broadcast_coords(ds, coords=("lon", "lat")):
    """broadcast 1d global coordinates"""
    lat2d, lon2d = xr.broadcast(ds[coords[1]], ds[coords[0]])
    lam, phi = lon2d, lat2d
    return lam, phi


def remap(ds, domain_info, vc, surflib, initial=False):
    """Remap CF-compliant GCM/analysis data to a REMO target domain.

    The workflow mirrors the legacy FORTRAN preprocessor (``intorg``):
    1. Horizontal interpolation of 3‑D and 2‑D fields from the input (global) grid
       to the rotated REMO grid using great‑circle search indices.
    2. First pressure adjustment using a planetary boundary layer (PBL) estimate.
    3. Vertical interpolation from input hybrid coefficients (``akgm``, ``bkgm``)
       to the requested target vertical coordinate ``vc`` (``ak``, ``bk``).
    4. Second pressure adjustment including hydrostatic consistency with geopotential.
    5. Wind rotation and divergence correction.
    6. Surface field remapping (SST, sea-ice) and optional initial-field enrichment.

    Parameters
    ----------
    ds : xarray.Dataset
        CF-compliant dataset holding at least the variables
        ``ta`` (air temperature), ``ua``/``va`` (winds), ``hus`` (specific humidity),
        ``ps`` (surface pressure), ``tos`` (sea surface temperature),
        ``orog`` (orography), ``sftlf`` (land fraction). Optionally ``clw`` for
        cloud liquid water content. Must also contain hybrid coordinate
        coefficients ``akgm`` and ``bkgm`` (full level) required for vertical interpolation.
    domain_info : dict
        Dictionary describing the target rotated domain. Required keys typically are
        ``ll_lon``, ``ll_lat`` (lower-left corner), ``dlon``, ``dlat`` (increment),
        ``nlon``, ``nlat`` (grid dimensions), ``pollon`` and ``pollat`` (rotated pole),
        plus optional metadata such as ``CORDEX_domain`` or ``domain_id``.
    vc : pandas.DataFrame
        Vertical coordinate definition with columns ``ak`` and ``bk`` (half levels)
        and ``akh`` / ``bkh`` (mid levels) as produced by ``pr.vc.tables``.
    surflib : xarray.Dataset
        Surface library on the target grid providing (at least) ``BLA`` (land/sea mask)
        and ``FIB`` (surface geopotential / orography divided by gravity). Additional
        static fields will be passed through if present.
    initial : bool, default: False
        If ``True`` add static/dynamic fields needed for initial condition creation
        (soil moisture/temperature etc.).

    Returns
    -------
    xarray.Dataset
        Forcing dataset on the target domain with rotated coordinates (``rlon``, ``rlat``)
        and vertical dimension ``lev``. Includes at least
        ``T`` (temperature), ``U`` / ``V`` (winds), ``PS`` (surface pressure),
        ``QD`` (specific humidity), ``QW`` (total water), ``QDBL`` (boundary layer
        humidity if produced), ``TSW`` (surface temperature over water), ``TSI`` (surface
        temperature over ice), ``SEAICE`` (fractional sea ice), vertical coordinate
        coefficients (``hyai``, ``hybi``, ``hyam``, ``hybm``) and auxiliary diagnostic
        attributes.

    Notes
    -----
    The function performs nearest-neighbour selection for final grid alignment after
    interpolation to ensure exact coordinate matching with ``surflib``. Calendar and
    time metadata are preserved from the input dataset where possible.
    """
    # rename vertical coordinate of input to avoid conflict with output lev
    gds = ds.copy()
    gds = gds.rename({gds.ta.cf["vertical"].name: lev_input})

    # remove time dimension if there is one
    surflib = surflib.squeeze(drop=True)
    surflib["rlon"] = np.round(surflib.rlon, 14)
    surflib["rlat"] = np.round(surflib.rlat, 14)

    fibem = surflib.FIB * const.grav_const
    blaem = surflib.BLA

    lamem, phiem = geo_coords(domain_info, surflib.rlon, surflib.rlat)

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
        clw = gds.clw
    else:
        clw = None

    arfgm = relative_humidity(gds.hus, gds.ta, gds.ps, gds.akgm, gds.bkgm, clw)
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
        blaem=blaem,
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
            blaem=blaem,
        )
    else:
        seaice = physics.seaice(tsw)

    water_content = physics.water_content(tem, arfem, psem, akbkem.akh, akbkem.bkh)
    tsi = physics.tsi(tsw)

    ads = xr.merge(
        [tem, uem_corr, vem_corr, psem, arfem, tsw, seaice, water_content, tsi, akbkem]
    )

    # rename for remo to recognize
    ads = ads.rename({"ak": "hyai", "bk": "hybi", "akh": "hyam", "bkh": "hybm"})

    # set global attributes
    ads.attrs = gds.attrs

    ads.attrs["history"] = "preprocessing with pyremo = {}".format(pr.__version__)
    ads.attrs["CORDEX_domain"] = domain_info.get("CORDEX_domain", "custom")
    ads.attrs["domain_id"] = domain_info.get("domain_id", "custom")

    # transpose to remo convention
    ads = ads.transpose(..., "lev", "rlat", "rlon")

    if initial is True:
        ads = add_surflib(ads, surflib)
        soil = _remap_era_soil(
            gds, tsw, fibem, blaem, lamem, phiem, lamgm, phigm, indii, indjj
        ).transpose("rlat", "rlon")
        # ads = xr.merge([ads, soil])
        ads = ads.merge(soil, join="override")
        ads["WS"] = np.minimum(ads.WSMX, ads.WS * ads.WSMX)
        ads["GLAC"] = xr.where(ads.SN > 9.5, 1.0, 0.0)

    grid = get_grid(domain_info)

    ads = ads.sel(rlon=grid.rlon, rlat=grid.rlat, method="nearest")
    ads["rlon"] = grid.rlon
    ads["rlat"] = grid.rlat
    ads = xr.merge([ads, grid])
    ads = update_attrs(ads)

    return ads


def remap_remo(
    ds,
    domain_info,
    vc,
    surflib,
    initial=False,
    lice=True,
    uvcor=True,
    domain_info_em=None,
):
    """Remap nested REMO (external model) output to a higher-resolution REMO domain.

    This is the double-nesting analogue of :func:`remap` operating on an *existing*
    REMO atmospheric forcing / restart (``t`` file) instead of CF GCM data. It reprojects
    the coarse external model (EM) REMO grid to a higher-resolution host model (HM) grid
    and optionally augments fields required for initial conditions.

    Parameters
    ----------
    ds : xarray.Dataset
        REMO external / driving model dataset containing at minimum 3‑D atmospheric
        fields (``T``, ``U``, ``V``, ``QD`` (specific humidity), ``QW`` (total water),
        ``PS``) and 2‑D surface fields (``TSW``, ``TSI``, ``FIB``, ``BLA``). May include
        ``SEAICE``; if absent and ``lice`` is True it will be diagnosed from ``TSW``.
    domain_info : dict
        Target (high-resolution) REMO domain description (see :func:`remap`).
    vc : pandas.DataFrame
        Target vertical coordinate definition with ``ak`` / ``bk`` (half) and
        ``akh`` / ``bkh`` (mid) coefficients.
    surflib : xarray.Dataset
        Surface library on the target domain providing at least ``BLA`` and ``FIB``.
    initial : bool, default: False
        If True, include soil fields and perform additional adjustments required
        for model cold / warm starts (e.g. soil moisture normalization, glacier mask).
    lice : bool or None, default: True
        If True, (re)compute sea ice fraction from surface temperature; if False, keep
        provided ``SEAICE`` values; if ``None`` auto-detect (derive when ``SEAICE`` missing).
    uvcor : bool, default: True
        Apply wind vector correction / rotation to remove artificial divergence and
        conform to the rotated pole geometry.
    domain_info_em : dict, optional
        Domain description of the external (coarser) REMO grid. If not supplied it is
        inferred from ``ds.cx.info()``.

    Returns
    -------
    xarray.Dataset
        Nested-domain forcing dataset on (``rlon``, ``rlat``) with vertical ``lev`` plus
        surface and (optionally) soil fields. Contains the transformed dynamic fields
        (``T``, ``U``, ``V``, ``PS``, humidity diagnostics, ``TSW``, ``TSI``, ``SEAICE``),
        vertical coordinate coefficients (``hyai``, ``hybi``, ``hyam``, ``hybm``) and
        any soil variables added for initialisation when ``initial`` is True.

    Notes
    -----
    The function performs two-stage pressure correction analogous to the global remap,
    but uses REMO-native hybrid coefficients from the external model (EM) as source.
    Wind fields are interpolated on staggered points and then re-centred / rotated.
    """
    tds = ds.copy()

    domain_info_hm = domain_info
    domain_info_em = domain_info_em or tds.cx.info()

    # check if we have to derive seaice from tsw
    has_seaice = "SEAICE" in tds
    if has_seaice is True:
        tds["SEAICE"] = tds.SEAICE.fillna(0.0)
    if lice is None:
        lice = "SEAICE" not in tds

    # rename coords so they dont conflict with hm coords
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

    indemi, indemj, dxemhm, dyemhm = intersect_regional(domain_info_em, domain_info_hm)
    indemi = indemi.assign_coords(rlon=surflib.rlon, rlat=surflib.rlat)
    indemj = indemj.assign_coords(rlon=surflib.rlon, rlat=surflib.rlat)
    dxemhm = dxemhm.assign_coords(rlon=surflib.rlon, rlat=surflib.rlat)
    dyemhm = dyemhm.assign_coords(rlon=surflib.rlon, rlat=surflib.rlat)

    # horizontal interpolation
    teh = interpolate_horizontal_remo(tds.T, indemi, indemj, dxemhm, dyemhm, "T")
    pseh = interpolate_horizontal_remo(tds.PS, indemi, indemj, dxemhm, dyemhm, "PS")
    pseh.name = "PSEH"
    ueh = interpolate_horizontal_remo(tds.U, indemi, indemj, dxemhm, dyemhm, "U", 1)
    veh = interpolate_horizontal_remo(tds.V, indemi, indemj, dxemhm, dyemhm, "V", 4)
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

    ficem = geopotential(tds.FIB, tds.T, tds.QD, tds.PS, tds.hyai, tds.hybi)

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

    # vertical interpolation
    akhem = tds.hyam
    bkhem = tds.hybm
    akbkhm = get_akbkem(vc)
    ptop = akbkhm.ak[0]

    thm = interpolate_vertical(
        teh, pseh, ps1hm, akhem, bkhem, akbkhm.akh, akbkhm.bkh, "T", kpbl, ptop=ptop
    )

    arfhm = interpolate_vertical(
        arfeh, pseh, ps1hm, akhem, bkhem, akbkhm.akh, akbkhm.bkh, "RF", kpbl, ptop=ptop
    )

    # second pressure correction and vertical interpolation of wind
    pshm = pressure_correction_ge(ps1hm, thm, arfhm, ficeh, fibhm, akbkhm.ak, akbkhm.bk)
    pshm.name = "PS"

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

    # set global attributes
    ads.attrs = tds.attrs

    ads.attrs["history"] = "preprocessing with pyremo = {}".format(pr.__version__)
    ads.attrs["CORDEX_domain"] = domain_info.get("CORDEX_domain", "custom")
    ads.attrs["domain_id"] = domain_info.get("domain_id", "custom")

    # transpose to remo convention
    return ads.transpose(..., "lev", "rlat", "rlon")


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
    # Reference: https://gitlab.dkrz.de/remo/RemapToRemo/-/blob/master/source/kernel/remo/bodfld.f90
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


def _remap_era_soil(ds, tswge, fibem, blaem, lamem, phiem, lamgm, phigm, indii, indjj):

    # based on https://gitlab.dkrz.de/remo/RemapToRemo/-/blob/master/setups/era_interim/bodfld.f90
    fak1 = 0.07
    fak2 = 0.21
    fak3 = 0.72
    # WSGm(ij) = (fak1*zwb1(ij)) + (fak2*zwb2(ij)) + (fak3*zwb3(ij))
    wsgm = fak1 * ds.swvl1 + fak2 * ds.swvl2 + fak3 * ds.swvl3

    # horizontal interpolation
    wsge = interpolate_horizontal(
        wsgm.clip(min=0.0), lamem, phiem, lamgm, phigm, "WS", indii=indii, indjj=indjj
    ).clip(min=0.0)
    td3ge = interpolate_horizontal(
        ds.tsl1, lamem, phiem, lamgm, phigm, "TD3", indii=indii, indjj=indjj
    )
    td4ge = interpolate_horizontal(
        ds.tsl2, lamem, phiem, lamgm, phigm, "TD4", indii=indii, indjj=indjj
    )
    td5ge = interpolate_horizontal(
        ds.tsl3, lamem, phiem, lamgm, phigm, "TD5", indii=indii, indjj=indjj
    )
    tdge = interpolate_horizontal(
        ds.tsl4, lamem, phiem, lamgm, phigm, "TD", indii=indii, indjj=indjj
    )
    snem = interpolate_horizontal(
        ds.snd, lamem, phiem, lamgm, phigm, "SN", indii=indii, indjj=indjj
    )
    wlem = interpolate_horizontal(
        ds.src.clip(min=0.0), lamem, phiem, lamgm, phigm, "WL", indii=indii, indjj=indjj
    ).clip(min=0.0)
    fibge = interpolate_horizontal(
        ds.orog, lamem, phiem, lamgm, phigm, "FIB", indii=indii, indjj=indjj
    )
    tslge = interpolate_horizontal(
        ds.skt, lamem, phiem, lamgm, phigm, "TSL", indii=indii, indjj=indjj
    )

    # wsmxem = surflib.WSMX.squeeze(drop=True)
    # wsem = np.minimum(wsmxem, wsge * wsmxem)

    # Initialize temperatures
    tslem, tswem, tsnem, td3em, td4em, td5em, tdem, tdclem = (
        physics.adapt_soil_temperatures(
            tdge, tswge, tslge, td3ge, td4ge, td5ge, fibem, fibge, blaem
        )
    )

    return xr.merge(
        [tslem, tswem, tsnem, td3em, td4em, td5em, tdem, tdclem, wlem, snem, wsge]
    ).squeeze(drop=True)


def remap_era_soil(ds, domain_info, surflib):
    """Create initial soil dataset from atmospheric dataset."""
    # rename vertical coordinate of input to avoid conflict with output lev
    ds = ds.copy()

    fak1 = 0.07
    fak2 = 0.21
    fak3 = 0.72
    # WSGm(ij) = (fak1*zwb1(ij)) + (fak2*zwb2(ij)) + (fak3*zwb3(ij))
    wsgm = fak1 * ds.swvl1 + fak2 * ds.swvl2 + fak3 * ds.swvl3

    # remove time dimension if there is one
    fibem = surflib.FIB.squeeze(drop=True) * const.grav_const
    blaem = surflib.BLA.squeeze(drop=True)

    lamem, phiem = geo_coords(domain_info, fibem.rlon, fibem.rlat)

    # broadcast 1d global coordinates
    lamgm, phigm = broadcast_coords(ds)

    print("getting matrix")
    # compute remap matrix
    indii, indjj = intersect(lamgm, phigm, lamem, phiem)  # .compute()

    print("remapping")
    # horizontal interpolation
    wsge = interpolate_horizontal(
        wsgm, lamem, phiem, lamgm, phigm, "WS", indii=indii, indjj=indjj
    )

    td3ge = interpolate_horizontal(
        ds.tsl1, lamem, phiem, lamgm, phigm, "TD3", indii=indii, indjj=indjj
    )
    td4ge = interpolate_horizontal(
        ds.tsl2, lamem, phiem, lamgm, phigm, "TD4", indii=indii, indjj=indjj
    )
    td5ge = interpolate_horizontal(
        ds.tsl3, lamem, phiem, lamgm, phigm, "TD5", indii=indii, indjj=indjj
    )
    tdge = interpolate_horizontal(
        ds.tsl4, lamem, phiem, lamgm, phigm, "TD", indii=indii, indjj=indjj
    )
    snem = interpolate_horizontal(
        ds.snd, lamem, phiem, lamgm, phigm, "SN", indii=indii, indjj=indjj
    )
    wlem = interpolate_horizontal(
        ds.src, lamem, phiem, lamgm, phigm, "WL", indii=indii, indjj=indjj
    )
    fibge = interpolate_horizontal(
        ds.orog, lamem, phiem, lamgm, phigm, "FIB", indii=indii, indjj=indjj
    )
    tslge = interpolate_horizontal(
        ds.skt, lamem, phiem, lamgm, phigm, "TSL", indii=indii, indjj=indjj
    )

    tswge = remap_sst(
        ds.tos,
        lamem,
        phiem,
        lamgm,
        phigm,
        blagm=xr.where(
            ds.tos.isel(time=0).isnull() if "time" in ds else ds.tos.isnull(), 1.0, 0.0
        ),
        blaem=blaem,
    )

    wsmxem = surflib.WSMX.squeeze(drop=True)
    wsem = np.minimum(wsmxem, wsge * wsmxem)

    # Initialize temperatures
    tslem, tswem, tsnem, td3em, td4em, td5em, tdem, tdclem = (
        physics.adapt_soil_temperatures(
            tdge, tswge, tslge, td3ge, td4ge, td5ge, fibem, fibge, blaem
        )
    )

    soil = xr.merge(
        [tslem, tswem, tsnem, td3em, td4em, td5em, tdem, tdclem, wlem, snem, wsem]
    )
    soil.attrs["history"] = "preprocessing with pyremo = {}".format(pr.__version__)
    soil.attrs["domain_id"] = domain_info.get("domain_id", "UNKNOWNs")

    soil = update_attrs(soil)

    # transpose to remo convention
    return soil.transpose(..., "rlat", "rlon")


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
    rlon_start = surflib.rlon.sel(rlon=ads.rlon.min(), method="nearest").item()
    rlon_end = surflib.rlon.sel(rlon=ads.rlon.max(), method="nearest").item()
    rlat_start = surflib.rlat.sel(rlat=ads.rlat.min(), method="nearest").item()
    rlat_end = surflib.rlat.sel(rlat=ads.rlat.max(), method="nearest").item()
    print(rlat_start, rlat_end, rlon_start, rlon_end)
    surflib = surflib.copy().sel(
        rlon=slice(rlon_start, rlon_end),
        rlat=slice(rlat_start, rlat_end),
    )
    # surflib = surflib.copy().isel(
    #     rlon=slice(1, -1),
    #     rlat=slice(1, -1),
    # )
    if "rotated_pole" in surflib:
        surflib = surflib.drop_vars("rotated_pole")
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
