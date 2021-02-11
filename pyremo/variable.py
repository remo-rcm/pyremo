#! /usr/bin/python
# coding: utf-8

"""This module is for handling variable meta information.


Example:
    Create an instance of :class:`RemoVariable`::

        from pyremo.variable import RemoVariable

        t = RemoVariable('T')
        print(t)
        print(t.__dict__)

"""

dump_str = """
 -------------
 Variable Info
 -------------
"""

import logging

import numpy as np
from .codes import get_dict
from . import physics
from .physics.units import ureg, Q_


class _Variable:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @property
    def properties(self):
        return self.__dict__.keys()

    def __str__(self):
        text = dump_str
        for prop in self.properties:
            text += " {:20}: {:20}\n".format(prop, str(getattr(self, prop)))
        return text


class RemoVariable(_Variable):
    """A Remo Variable meta information container.

    The attributes of an instance of RemoVariable dependent on the
    code list content. The attributes will be determined by the
    entries found for the variable idenfitifer.

    Attributes:
        variable (str): Variable name.
        code (int): Code.

    """

    def __init__(self, id, data=None, **kwargs):
        """Creates an instance.

        Args:
            id (int or str): Variable identifier.
        """
        _Variable.__init__(self, **kwargs)
        self.id = id
        self._get_info(id)
        self.data = self._set_data(data)

    def _get_info(self, id):
        self.varinfo = get_dict(id)
        if self.varinfo:
            for prop, value in self.varinfo.items():
                setattr(self, prop, value)
        else:
            raise Exception("Undefined Identifier: {}!".format(id))

    def _set_data(self, data):
        if data is not None and not isinstance(data, VariableData):
            if not hasattr(data, "units"):
                return VariableData(data, self.units)
            else:
                return VariableData(data)
        return data

    def to(self, units):
        """convert units of data."""
        data = self.data.to(units)
        remo_var = RemoVariable(self.id, data)
        remo_var.units = str(units)
        return remo_var


class VariableData:
    """A Data container for Variables."""

    def __init__(self, data, units=None, **kwargs):
        self.data = data
        if units is None and hasattr(self.data, "units"):
            self._set_units(self.data.units)
        else:
            self._set_units(units)

    def __str__(self):
        return str(self.data)

    def _set_units(self, units):
        if isinstance(units, type(ureg.unit)):
            self.units = units
        elif units:
            self.units = physics.units.units(units)

    def to(self, units):
        """convert units of data."""
        data = self.data.copy()
        data[:] = Q_(np.array(self.data), self.units).to(units)
        if hasattr(data, "units"):
            data.attrs["units"] = str(units)
        return VariableData(data, units)


def to(data, units, in_units=None):
    """convert units of data"""
    if in_units is None and hasattr(data, "units"):
        in_units = data.units
    if in_units:
        data = data.copy()
        if hasattr(data, "units"):
            data.attrs["units"] = str(units)
        data[:] = Q_(np.array(data), in_units).to(units)
        return data
    raise Exception("could not determine unit of data")


def remo_variable(data, id):
    pass


def variable(id):
    """Top level function to create a Remo variable.

    Args:
        id (int or str): Variable identifier.

    """
    return RemoVariable(id)
