"""Top-level package for pyremo."""

__author__ = """Lars Buntemeyer"""
__email__ = "lars.buntemeyer@hzg.de"


# from . import core, preproc
from . import preproc, physics

from .core.domain import domain_info, remo_domain
from .core.remo_ds import open_remo_dataset, update_meta_info, preprocess
from .core.cal import parse_dates, parse_absolute_time
from .core import remo_tutorial as tutorial, remo_ds, data
from .tables import domains, vc
from .core import codes
from .cmor.cmor import cmorize_variable


from .version import version

__version__ = version
