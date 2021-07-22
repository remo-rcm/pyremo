"""Top-level package for pyremo."""

__author__ = """Lars Buntemeyer"""
__email__ = "lars.buntemeyer@hzg.de"
__version__ = "0.1.0"


from . import core, preproc

from .core.domain import domain_info, remo_domain
from .core.remo_ds import open_remo_dataset, parse_dates, update_meta_info
from .core import data
from .tables import domains, codes, vc

