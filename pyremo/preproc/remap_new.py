import cordex as cx
import xarray as xr
from xesmf import Regridder

import pyremo as pr

from . import physics
from .constants import lev_input
from .core import (
    correct_uv,
    geopotential,
    get_akbkem,
    interpolate_vertical,
    pbl_index,
    pressure_correction_em,
    pressure_correction_ge,
    relative_humidity,
    rotate_uv,
)


class const:
    """constants used for unit conversion"""

    grav_const = 9.806805923
    absolute_zero = 273.5


def get_grid(domain_info):
    return cx.create_dataset(**domain_info)


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
        except Exception:
            pass
    return ds


def remap(
    gds,
    domain_info,
    vc,
    surflib,
    method="bilinear",
    periodic=True,
    weights=None,
    weights_u=None,
    weights_v=None,
    weights_with_mask=None,
):
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

    #    domain_info["ll_lon"] = domain_info["ll_lon"] - 1.0 * domain_info["dlon"]
    #    domain_info["ll_lat"] = domain_info["ll_lat"] - 1.0 * domain_info["dlat"]
    #    domain_info["nlon"] = domain_info["nlon"] + 2
    #    domain_info["nlat"] = domain_info["nlat"] + 2
    # regional grid
    grid = get_grid(domain_info)

    domain_info_u = domain_info.copy()
    domain_info_v = domain_info.copy()
    domain_info_u["ll_lon"] = domain_info["ll_lon"] + 0.5 * domain_info["dlon"]
    domain_info_v["ll_lat"] = domain_info["ll_lat"] + 0.5 * domain_info["dlon"]
    grid_u = get_grid(domain_info_u)
    grid_v = get_grid(domain_info_v)
    # return grid_u, grid_v
    # remove time dimension if there is one
    surflib = surflib.copy().isel(rlon=slice(1, -1), rlat=slice(1, -1))
    fibem = (
        surflib.FIB.squeeze(drop=True).assign_coords(rlon=grid.rlon, rlat=grid.rlat)
        * const.grav_const
    )

    # lamem, phiem = geo_coords(domain_info, fibem.rlon, fibem.rlat)
    lamem = xr.concat([grid.lon, grid_u.lon, grid_v.lon], dim="pos", join="override")
    phiem = xr.concat([grid.lat, grid_u.lat, grid_v.lat], dim="pos", join="override")
    # lamem, phiem = geo_coords(domain_info, fibem.rlon, fibem.rlat)#.isel(rlon=slice(1,-1), rlat=slice(1,-1))
    # grid_u = grid_u.rename({"rlon":"rlon_u"})
    # grid_v = grid_v.rename({"rlat":"rlat_v"})

    # return lamem
    # broadcast 1d global coordinates
    # lamgm, phigm = broadcast_coords(gds)

    # compute remap matrix
    # indii, indjj = intersect(lamgm, phigm, lamem, phiem)  # .compute()

    regridder = Regridder(
        gds,
        grid,
        method=method,
        periodic=periodic,
        reuse_weights=(weights is not None),
        weights=weights,
    )
    regridder_u = Regridder(
        gds,
        grid_u,
        method=method,
        periodic=periodic,
        reuse_weights=(weights_u is not None),
        weights=weights_u,
    )
    regridder_v = Regridder(
        gds,
        grid_v,
        method=method,
        periodic=periodic,
        reuse_weights=(weights_v is not None),
        weights=weights_v,
    )
    print(surflib.BLA)
    regridder_sst = Regridder(
        ~gds.tos.isnull().squeeze(drop=True).to_dataset().rename({"tos": "mask"}),
        xr.merge([grid, 1.0 - surflib.BLA.squeeze(drop=True)], join="override").rename(
            {"BLA": "mask"}
        ),
        method=method,
        periodic=periodic,
        reuse_weights=(weights_with_mask is not None),
        weights=weights_with_mask,
        extrap_method="inverse_dist",
        ignore_degenerate=True,
    )
    print(regridder_sst)

    tge = regridder(gds.ta)
    psge = regridder(gds.ps)
    uge = regridder(gds.ua)
    vge = regridder(gds.ua)
    uvge = regridder_v(gds.ua)
    vuge = regridder_u(gds.va)
    fibge = regridder(gds.orog)

    uvge = uvge.assign_coords(rlon=tge.rlon, rlat=tge.rlat)
    vuge = vuge.assign_coords(rlon=tge.rlon, rlat=tge.rlat)

    ficgm = geopotential(
        gds.orog, gds.ta, gds.hus, gds.ps, gds.akgm, gds.bkgm
    )  # .squeeze(drop=True)

    ficge = regridder(ficgm)
    # ficge = interpolate_horizontal(
    #     ficgm, lamem, phiem, lamgm, phigm, "FIC", indii=indii, indjj=indjj
    # )

    if "clw" in gds:
        # if False:
        arfgm = relative_humidity(gds.hus, gds.ta, gds.ps, gds.akgm, gds.bkgm, gds.clw)
    else:
        arfgm = relative_humidity(gds.hus, gds.ta, gds.ps, gds.akgm, gds.bkgm)
    arfge = regridder(arfgm)
    # arfge = interpolate_horizontal(
    #     arfgm, lamem, phiem, lamgm, phigm, "AREL HUM", indii=indii, indjj=indjj
    # )

    # return uge, vge, uvge, vuge, lamem, phiem
    # wind vector rotation
    uge_rot, vge_rot = rotate_uv(
        uge, vge, uvge, vuge, lamem, phiem, domain_info["pollon"], domain_info["pollat"]
    )
    # first pressure correction
    kpbl = pbl_index(gds.akgm, gds.bkgm)
    # return psge, tge, arfge, fibge, fibem
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
    philuem = domain_info["ll_lat"]  # + 1.0 * dphiem
    uem_corr, vem_corr = correct_uv(
        uem, vem, psem, akbkem.ak, akbkem.bk, lamem, phiem, philuem, dlamem, dphiem
    )

    tsw = regridder_sst(gds.tos)
    tsw.name = "TSW"

    # check if gcm contains seaice, else derive from sst
    if "sic" in gds:
        seaice = regridder_sst(gds.sic)
    else:
        seaice = physics.seaice(tsw)
    seaice.name = "SEAICE"

    water_content = physics.water_content(tem, arfem, psem, akbkem.akh, akbkem.bkh)
    tsi = physics.tsi(tsw)

    ads = xr.merge(
        [tem, uem_corr, vem_corr, psem, arfem, tsw, seaice, water_content, tsi, akbkem]
    )

    grid = get_grid(domain_info)

    # ads = ads.sel(rlon=grid.rlon, rlat=grid.rlat, method="nearest")
    # ads = ads.isel(rlon=slice(1,-1), rlat=slice(1,-1))
    ads["rlon"] = grid.rlon  # .isel(rlon=slice(1,-1))
    ads["rlat"] = grid.rlat  # .isel(rlat=slice(1,-1))

    ads = xr.merge([ads, grid])

    # rename for remo to recognize
    ads = ads.rename({"ak": "hyai", "bk": "hybi", "akh": "hyam", "bkh": "hybm"})

    # set global attributes
    ads.attrs = gds.attrs

    ads.attrs["history"] = "preprocessing with pyremo = {}".format(pr.__version__)
    ads.attrs["CORDEX_domain"] = domain_info.get("short_name", "no name")

    ads = update_attrs(ads)

    # transpose to remo convention
    return ads.transpose(..., "lev", "rlat", "rlon")


class Remapper:
    def __init__(
        self,
        gds,
        domain_info,
        vc,
        surflib,
        method_atmo="bilinear",
        # method_sst="nearest_s2d",
        method_sst="bilinear",
        periodic=True,
        sst=None,
    ):
        self.domain_info = domain_info
        self.vc = vc
        self.surflib = surflib.load()
        self.method = method_atmo
        self.method_sst = method_sst
        self.periodic = periodic
        self._create_rcm_grids(domain_info)
        # if "time" in gds.coords:
        self._init_weights(gds, sst=sst)
        # else:
        # self._init_regridder(gds.isel(time=0))

    def _create_rcm_grids(self, domain_info):
        # regional grid
        self.grid = get_grid(domain_info)
        domain_info_u = domain_info.copy()
        domain_info_v = domain_info.copy()
        domain_info_u["ll_lon"] = domain_info["ll_lon"] + 0.5 * domain_info["dlon"]
        domain_info_v["ll_lat"] = domain_info["ll_lat"] + 0.5 * domain_info["dlon"]
        self.grid_u = get_grid(domain_info_u)
        self.grid_v = get_grid(domain_info_v)

    def _init_weights(self, gds, sst=None):
        self.weights = Regridder(
            gds, self.grid, method=self.method, periodic=self.periodic
        ).to_netcdf("weights.nc")
        self.weights_u = Regridder(
            gds, self.grid_u, method=self.method, periodic=self.periodic
        ).to_netcdf("weights_u.nc")
        self.weights_v = Regridder(
            gds, self.grid_v, method=self.method, periodic=self.periodic
        ).to_netcdf("weights_v.nc")
        if sst is None:
            sst = gds.tos
        sst_mask = ~sst.isnull().squeeze(drop=True)
        sst_mask.name = "mask"
        bla = self.surflib.BLA.isel(rlon=slice(1, -1), rlat=slice(1, -1))
        bla.name = "mask"
        self.weights_sst = Regridder(
            sst_mask.to_dataset(),
            xr.merge([self.grid, 1.0 - bla], join="override"),
            # self.grid.merge(1 - bla, join="override"),
            method=self.method_sst,
            extrap_method="inverse_dist",
            ignore_degenerate=True,
            periodic=self.periodic,
        ).to_netcdf("weights_sst.nc")

    def remap(self, gds):
        return remap(
            gds,
            self.domain_info,
            self.vc,
            self.surflib,
            weights=self.weights,
            weights_u=self.weights_u,
            weights_v=self.weights_v,
            weights_with_mask=self.weights_sst,
        )
