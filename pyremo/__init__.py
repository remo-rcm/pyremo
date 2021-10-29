"""Top-level package for pyremo."""

__author__ = """Lars Buntemeyer"""
__email__ = "lars.buntemeyer@hzg.de"


# from . import core, preproc
from . import preproc, physics

from .core.domain import domain_info, remo_domain
from .core.remo_ds import open_remo_dataset, parse_dates, update_meta_info
from .core import remo_tutorial as tutorial, remo_ds, data
from .tables import domains, vc
from .core import codes


from .version import version

__version__ = version
