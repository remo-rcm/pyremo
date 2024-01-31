"""
Useful for:
* users learning pyremo
* building tutorials in the documentation.
* stolen from xarray!
"""

# code stolen from xarray, I am sorry!
import os
import pathlib

import xarray as xr
from xarray import open_dataset as _open_dataset

base_url = "https://github.com/remo-rcm/pyremo-data"
version = "main"

_default_cache_dir_name = "pyremo_tutorial_data"

external_urls = {}  # type: dict

file_formats = {}


def _construct_cache_dir(path):
    import pooch

    if isinstance(path, os.PathLike):
        path = os.fspath(path)
    elif path is None:
        path = pooch.os_cache(_default_cache_dir_name)

    return path


def _check_netcdf_engine_installed(name):
    version = file_formats.get(name)
    if version == 3:
        try:
            import scipy  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either scipy or "
                    "netCDF4 to be installed."
                )
    if version == 4:
        try:
            import h5netcdf  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either h5netcdf "
                    "or netCDF4 to be installed."
                )


# idea borrowed from Seaborn
def open_dataset(
    name="remo_EUR-44",
    cache=True,
    cache_dir=None,
    *,
    engine=None,
    **kws,
):
    """
    Open a dataset from the online repository (requires internet).

    If a local copy is found then always use that to avoid network traffic.

    Available datasets:

    * ``"remo_EUR-44"``: remo example output on EUR-44 domain
    * ``"remo_EUR-11_TEMP2"``: remo 2m temperature time series on EUR-11 domain

    Parameters
    ----------
    name : str
        Name of the file containing the dataset.
        e.g. 'remo_EUR-44'
    cache_dir : path-like, optional
        The directory in which to search for and write cached data.
    cache : bool, optional
        If True, then cache data locally for use on subsequent calls
    **kws : dict, optional
        Passed to xarray.open_dataset

    See Also
    --------
    xarray.open_dataset
    """
    try:
        import pooch
    except ImportError as e:
        raise ImportError(
            "tutorial.open_dataset depends on pooch to download and manage datasets."
            " To proceed please install pooch."
        ) from e

    logger = pooch.get_logger()
    logger.setLevel("WARNING")

    cache_dir = _construct_cache_dir(cache_dir)
    if name in external_urls:
        url = external_urls[name]
    else:
        path = pathlib.Path(name)
        if not path.suffix:
            # process the name
            default_extension = ".nc"
            if engine is None:
                _check_netcdf_engine_installed(name)
            path = path.with_suffix(default_extension)
        elif path.suffix == ".grib":
            if engine is None:
                engine = "cfgrib"

        url = f"{base_url}/raw/{version}/{path.name}"

    # retrieve the file
    filepath = pooch.retrieve(url=url, known_hash=None, path=cache_dir)
    ds = _open_dataset(filepath, engine=engine, **kws)
    if not cache:
        ds = ds.load()
        pathlib.Path(filepath).unlink()

    return ds


def load_dataset(*args, **kwargs):
    """
    Open, load into memory, and close a dataset from the online repository
    (requires internet).

    See Also
    --------
    open_dataset
    """
    with open_dataset(*args, **kwargs) as ds:
        return ds.load()


def mpi_esm(**kwargs):
    files = [
        "ta_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "ua_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "va_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "hus_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "ps_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "orog_fx_MPI-ESM1-2-HR_historical_r1i1p1f1_gn",
        "sftlf_fx_MPI-ESM1-2-HR_historical_r1i1p1f1_gn",
    ]

    return xr.merge(
        [load_dataset(f, **kwargs) for f in files], compat="override", join="override"
    )


def mpi_esm_tos(**kwargs):
    return open_dataset(
        "tos_Oday_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_19781231-19790102", **kwargs
    )
