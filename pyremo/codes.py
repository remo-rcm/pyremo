#! /usr/bin/python
# coding: utf-8
# flake8: noqa

"""This module is for defining the Remo Code Table

The table in this module contain some meta information about REMO codes
and variable names. The tables define variable names
and codes for the REMO2015 version. The mostly contains:

.. data:: table

**code:**
    integer code.
**name:**
    variable name (REMO standard)
**cf_name:**
    variable name fullfilling Climate and Forecast convention (if available)
**long_name:**
    description (REMO standard)
**unit:**
    unit (REMO standard)
**time_cell_method:**
    time cell method for standard Remo output (mean or inst)


Example:
    Work with the code table and retrieve some information ,e.g.,::

        from pyremo.codes import table, get_dict

        # print the whole table
        print(table)
        # get dicitonary of a single variable
        t = get_dict('T')
        print(t)
        # or use the code...
        t = get_dict(130)


"""

# from .tables import codes as code_table  # code_table, read_table

import numpy as np
import pandas as pd

from .tables import codes

# table = pd.concat([table for name, table in code_table.items()])


def get_dict(id):
    """Returns a dictionary with variable info.

    Searches the code table for a certain variable id.

    Args:
        id (int or str): The variable identifier (might be code or
            variable name or a variable name containing a code or
            the variable name according to CF conventions).

    Returns:
        varinfo (dict): Dictionary with variable information.

    """
    # id is expected to be a variable name
    if isinstance(id, (str)):
        try:
            return get_dict_by_name(id)
        except:
            try:
                return get_dict_by_code(_code_from_varname(id))
            except:
                raise Exception("unknown identifier: {}".format(id))
    # id is expected to be a code number
    elif np.issubdtype(type(id), np.integer):
        return get_dict_by_code(id)
    else:
        return None


def get_dict_by_name(varname):
    """Returns a dictionary with variable info.

    Searches the code table for a certain variable name.

    Args:
        varname (str): The variable name. If the name is not found
        in the table, it searches for the cf name.

    Returns:
        varinfo (dict): Dictionary with variable information.

    """
    table = pd.concat(codes.tables.values())
    df = table.loc[table["variable"] == varname]
    if df.empty:
        # try with cf name
        df = table.loc[table["cf_name"] == varname]
    code = df.index[0]
    dict = df.where(pd.notnull(df), None).to_dict(orient="list")
    dict = {key: item[0] for key, item in dict.items()}
    dict["code"] = code
    return dict


def get_dict_by_code(code):
    """Returns a dictionary with variable info.

    Searches the code table for a certain variable code.

    Args:
        code (int): The variable code.

    Returns:
        varinfo (dict): Dictionary with variable information.

    """
    table = pd.concat(codes.tables.values())
    series = table.loc[code]
    dict = series.where(pd.notnull(series), None).to_dict()
    dict["code"] = code
    return dict


def _code_from_varname(varname):
    """Returns a code from a varname.

    Used typically for varnames create by a `cdo -f nc copy` command
    on IEG files.
    """
    if "var" in varname:
        import re

        return int(re.findall("[0-9]+", varname)[0])
    else:
        return None
