import pkg_resources

from . import physics, tutorial
from .core import codes, data, remo_ds
from .core.cal import parse_absolute_time, parse_dates
from .core.conventions import output_pattern
from .core.domain import domain_info, magic_number, magic_numbers, remo_domain
from .core.remo_ds import open_remo_dataset, preprocess, update_meta_info
from .tables import domains, vc

try:
    __version__ = pkg_resources.get_distribution("pyremo").version
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"

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
    "magic_number",
    "magic_numbers",
    "update_meta_info",
    "domains",
    "vc",
    "output_pattern",
    "",
]
