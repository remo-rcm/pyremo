#! /usr/bin/python
# coding: utf-8
#
# This file is part of PyRemo. PyRemo is a toolbox to facilitate
# treatment and plotting of REMO or other rotated and non-rotated
# data.
#
# Copyright (C) 2010-2014 REMO Group
# See COPYING file for copying and redistribution conditions.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


"""This module is for defining the Remo Code Table


.. note::
    This module is autmatically created from a csv table using.


The table in this module contain some meta information about REMO codes
and variable names. There is only one table right now defining variable names
and codes for the REMO2015 version. The table contains right now:


.. data:: REMO2015

**name:**
    variable name (REMO standard)
**cf_name:**
    variable name fullfilling Climate and Forecast convention (if available)
**long_name:**
    description (REMO standard)
**unit:**
    unit (REMO standard)


"""

from ._tables import tables#code_table, read_table

import pandas as pd

table = pd.concat([table for name,table in tables['code'].items()])


def get_dict(id):
    """returns a dictionary with variable info.

    nan values will be replaced with None.
    """
    # id is expected to be a variable name
    if isinstance(id,(str)):
        return get_dict_by_name(id)
    # id is expected to be a code number
    elif isinstance(id, int):
        return get_dict_by_code(id)
    else:
        return None


def get_dict_by_name(varname):
    """returns a dictionary with variable info.

    nan values will be replaced with None.
    """
    df = table.loc[table['variable']==varname]
    code = df.index[0]
    dict = df.where(pd.notnull(df), None).to_dict(orient='list')
    dict = { key:item[0] for key, item in dict.items() }
    dict['code'] = code
    return dict


def get_dict_by_code(code):
    """returns a dictionary with variable info.

    nan values will be replaced with None.
    """
    series = table.loc[code]
    dict = series.where(pd.notnull(series), None).to_dict()
    dict['code'] = code
    return dict
