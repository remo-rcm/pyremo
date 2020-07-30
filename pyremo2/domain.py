#! /usr/bin/python
# coding: utf-8
#
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


from itertools import chain, product
from cordex.domain import Domain, domain_from_table


from . import _tables


tables = _tables.tables["domain"]


class _DomainFactory(object):
    """Factory class for creating a domain instance."""

    @classmethod
    def names_from_csv(cls, table=None):
        """Returns a list of names of csv domains."""
        if table:
            return list(tables[table].index.values)
        else:
            return list(
                chain.from_iterable(
                    [[sn for sn in t.index.values] for n, t in tables.items()]
                )
            )

    @classmethod
    def create_domain_from_table(cls, short_name):
        """Returns a list of names of csv domains."""
        for table_name, table in tables.items():
            if short_name in table.index.values:
                return domain_from_table(short_name, table)

    @classmethod
    def names(cls, table=None):
        """Returns a list of names of available Domains."""
        return cls.names_from_csv(table)

    @classmethod
    def domains(cls, table=None):
        """Returns a dictionary of names and domains."""
        return {name: cls.get_domain(name) for name in cls.names_from_csv(table)}

    @classmethod
    def get_domain(cls, short_name):
        """Returns a Domain instance.

        Args:
          name (str): standard name of the Domain.

        Returns:
          Domain (:class:`Domain`) : a Domain instance.

        """
        out = None
        if short_name in cls.names_from_csv():
            out = cls.create_domain_from_table(short_name)
        if out is None:
            _logger.error("Unknown domain name: " + short_name)
            _logger.info("Known domain names: " + str(cls.names()))
            raise Exception("Unknown domain name: " + short_name)
        else:
            return out


def domain(name):
    """Top level Domain function to get a :class:`Domain` instance.

    Args:
      name (str): name of the domain instance.

    Returns:
      :class:`Domain`: preconfigured domain instance.

    """
    return _DomainFactory().get_domain(name)


def domains(table=None):
    """Top level function that returns a dictionay of CORDEX domains.

    Returns:
      domains (dict): dict of available CORDEX domain names.

    """
    return _DomainFactory().domains(table)


def names(table=None):
    """Top level function that returns the names of available CORDEX domains.

    Returns:
      names (list): list of available CORDEX domain names.

    """
    return _DomainFactory().names(table)


# def table(name):
#    """Top level function that returns a CORDEX table.
#
#    Args:
#      name (str): name of the CORDEX table.
#
#    Returns:
#      table (DataFrame): Cordex table.
#
#    """
#    return table[name]

# def tables():
#    """Top level function that returns a list of available CORDEX tables.
#
#    Returns:
#      names (list): list of available CORDEX domains.
#
#    """
#    return list(table.keys())


def magic_number(n=1, m=0, o=0):
    """returns a magic number for REMO grid boxes."""
    return 2 ** n * 3 ** m * 5 ** o + 1


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
