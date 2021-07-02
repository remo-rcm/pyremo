# flake8: noqa
"""This module defines the csv tables for remo.
"""

import pandas as pd
import pkg_resources

from . import domains as dm
from . import code_list as cd
from . import vc as vcoord


def read_resource_table(resource, csv, **kwargs):
    """reads a csv table from the package resource."""
    csv_file = pkg_resources.resource_stream(resource, csv)
    return pd.read_csv(csv_file, **kwargs)


def read_resource_tables(resource, csv_dict, **kwargs):
    """reads all csv tables from the package resource."""
    tables = {}
    for set_name, table in csv_dict.items():
        tables[set_name] = read_resource_table(resource, table, **kwargs)
    return tables


domains = read_resource_tables(
    "pyremo.tables.domains", dm.tables, index_col="short_name"
)
codes = read_resource_tables("pyremo.core.tables.code_list", cd.tables, index_col="code")
vc = read_resource_tables("pyremo.tables.vc", vcoord.tables)
