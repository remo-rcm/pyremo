# flake8: noqa

from .cf import get_gcm_gfile, get_gcm_dataset, to_cfdatetime, gfile, get_gfile

# from .core import gfile
from .era5 import ERA5
from .preprocessor import (
    ERA5Preprocessor,
    CFPreprocessor,
    RemoPreprocessor,
)
from .remapping import remap, remap_remo, to_netcdf, to_tar

# from .remap_new import Remapper
from .utils import write_forcing_file


PREPROCESSOR_CLASSES = {
    "CF": CFPreprocessor,
    "ERA5": ERA5Preprocessor,
    "REMO": RemoPreprocessor,
    #  "cloud": CloudPreprocessor,
}


def get_preprocessor(preprocessor_name):
    """
    Get the preprocessor class based on the name.

    Args:
        preprocessor_name (str): Name of the preprocessor.

    Returns:
        class: Preprocessor class.

    Raises:
        ValueError: If the preprocessor name is not recognized.
    """
    try:
        return PREPROCESSOR_CLASSES[preprocessor_name]
    except KeyError:
        raise ValueError(f"Preprocessor '{preprocessor_name}' not recognized.")
