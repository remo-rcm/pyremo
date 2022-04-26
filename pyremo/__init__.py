"""Top-level package for pyremo."""

__author__ = """Lars Buntemeyer"""
__email__ = "lars.buntemeyer@hzg.de"


# from . import core, preproc
# from . import preproc, physics, cmor
from . import physics, tutorial
from .core import codes, data, remo_ds
from .core.cal import parse_absolute_time, parse_dates
from .core.domain import domain_info, remo_domain
from .core.remo_ds import open_remo_dataset, preprocess, update_meta_info
from .tables import domains, vc
from .version import version

__version__ = version

__all__ = [
    "physics",
    "tutorial",
    "codes",
    "data",
    "remo_ds",
    "parse_absolute_time",
    "parse_dates",
    "domain_info",
    "remo_domain",
    "open_remo_dataset",
    "preprocess",
    "update_meta_info",
    "domains",
    "vc",
    "version",
    "",
]
