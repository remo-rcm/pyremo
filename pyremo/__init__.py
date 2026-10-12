from importlib.metadata import version as _get_version

from . import codes, data, physics, remo_ds, tutorial
from .cal import parse_absolute_time, parse_dates
from .conventions import file_pattern, output_pattern
from .domain import domain_info, magic_number, magic_numbers, remo_domain
from .remo_ds import open_remo_dataset, preprocess, update_meta_info
from .tables import domains, vc

try:
    __version__ = _get_version("pyremo")
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"

__all__ = [
    "",
    "codes",
    "data",
    "domain_info",
    "domains",
    "file_pattern",
    "magic_number",
    "magic_numbers",
    "open_remo_dataset",
    "output_pattern",
    "parse_absolute_time",
    "parse_dates",
    "physics",
    "preprocess",
    "remo_domain",
    "remo_ds",
    "tutorial",
    "update_meta_info",
    "vc",
]
