#! /usr/bin/python
# coding: utf-8
# flake8: noqa
# This file is part of PyRemo. PyRemo is a toolbox to facilitate
# treatment and plotting of REMO or other rotated and non-rotated
# data.
#
# Copyright (C) 2010-2020 REMO Group
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
"""Domain module

This module creates and defines cordex domains for remo. It is based on the
cordex domain module.


.. note::
       A remo domain is usually a little larger than the "official" cordex
       domain according to the archive specification since the regional model
       usually accounts for "buffer" zone where the lateral boundary conditions
       are nudged.

"""

from itertools import chain, product

import cordex as cx
import pandas as pd

from .tables import domains


def remo_domain(domain_id, **kwargs):
    """Creates an xarray dataset containg the domain grid definitions.

    Parameters
    ----------
    domain_id:
        Name of the Cordex Domain.
    dummy : str or logical
        Name of dummy field, if dummy=topo, the cdo topo operator will be
        used to create some dummy topography data. dummy data is useful for
        looking at the domain with ncview.

    Returns
    -------
    Dataset : xarray.core.Dataset
        Dataset containing the coordinates.

    """
    return cx.cordex_domain(domain_id, tables=domains.table, **kwargs)


def domain_info(domain_id):
    """Returns a dictionary containg the domain grid definitions.

    Returns a dictionary with grid information according to the
    Cordex archive specifications.

    Parameters
    ----------
    domain_id:
        Name of the Cordex Domain.

    Returns
    -------
    domain info : dict
        Dictionary containing the grid information.

    """
    return cx.domain_info(domain_id, tables=domains.table)


def magic_number(n=1, m=0, o=0):
    """returns a magic number for REMO grid boxes."""
    return 2**n * 3**m * 5**o + 1


def magic_numbers(size=100):
    """returns a list magic numbers

    Args:
      size (int): length of list of magic numbers.

    Returns:
      magic numbers (list): list of magic numbers.
    """
    n = range(1, 15)
    m = range(0, 10)
    o = range(0, 10)
    numbers = [magic_number(*x) for x in product(n, m, o)]
    numbers.sort()
    return numbers[0:size]
