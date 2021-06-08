
import os

DKRZ_URL = "https://swift.dkrz.de/v1/dkrz_ffd3ca9004324ad28243244b834f92b1/remo/data"

flake_EUR_11_glcc = os.path.join(DKRZ_URL, 'flake/flake_v3_glcc_defD10.0m_frac_EUR-11.nc')
flake_EUR_44_glcc = os.path.join(DKRZ_URL, 'flake/flake_v3_glcc_defD10.0m_frac_EUR-44.nc')
tutorial_data = os.path.join(DKRZ_URL, 'example/e056111t2006010100.nc')

bodlib_tpl = os.path.join(DKRZ_URL, 'surface-library/lib_{}_frac')


from .remo_ds import open_remo_dataset


def get_file(url):
    import fsspec
    import xarray as xr
    with fsspec.open(url) as f:
        ds = xr.open_dataset(f)
    return ds


def bodlib(domain='EUR-11'):
    import fsspec
    url = bodlib_tpl.format(domain)
    with fsspec.open(url) as f:
        ds = open_remo_dataset(f)
    return ds


def example_eur44():
    import fsspec
    url = bodlib_tpl.format(domain)
    with fsspec.open(url) as f:
        ds = open_remo_dataset(f)
    return ds


def example_output():
    from . import remo_ds as rds
    return rds.update_meta_info(get_file(tutorial_data))
