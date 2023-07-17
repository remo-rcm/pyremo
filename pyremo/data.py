import os

from . import domain as dm
from . import remo_ds as rds

DKRZ_URL = "https://swift.dkrz.de/v1/dkrz_ffd3ca9004324ad28243244b834f92b1/remo/data"

flake_EUR_11_glcc = os.path.join(
    DKRZ_URL, "flake/flake_v3_glcc_defD10.0m_frac_EUR-11.nc"
)
flake_EUR_44_glcc = os.path.join(
    DKRZ_URL, "flake/flake_v3_glcc_defD10.0m_frac_EUR-44.nc"
)
tutorial_data = os.path.join(DKRZ_URL, "example/e056111t2006010100.nc")

bodlib_tpl = os.path.join(DKRZ_URL, "surface-library/lib_{}_frac.nc")


def _get_file(url):
    import fsspec
    import xarray as xr

    with fsspec.open(url) as f:
        ds = xr.open_dataset(f)
    return ds


def surflib(domain, crop=True, update_meta=True):
    """Return a REMO surface library from the DKRZ cloud.

    Please note that by default the surface library will be cropped to
    the size of the REMO domain. The original libraries contained
    an original boundary row.

    Requires network connection for online access.

    Parameters
    ----------
    domain : str
        REMO domain short name.
    crop : bool
        If True, crop to REMO domain.
    update_meta : bool
        If True, update meta information if required.

    Returns
    -------
    dataset : xr.Dataset
        Surface library dataset.

    """
    import fsspec
    import xarray as xr

    url = bodlib_tpl.format(domain)
    with fsspec.open(url) as f:
        # ds = open_remo_dataset(f, update_meta=True).squeeze(drop=True)
        ds = xr.open_dataset(f).squeeze(drop=True)
    if update_meta:
        ds = rds.update_meta_info(ds)
    if crop:
        grid = dm.remo_domain(domain)
        ds = ds.sel(rlon=grid.rlon, rlat=grid.rlat, method="nearest")
    return ds


# def example_eur44():
#    import fsspec
#    url = bodlib_tpl.format(domain)
#    with fsspec.open(url) as f:
#        ds = open_remo_dataset(f)
#    return ds


def example_output():
    """Returns a dataset containing REMO example output."""
    from . import remo_ds as rds

    return rds.update_meta_info(_get_file(tutorial_data))
