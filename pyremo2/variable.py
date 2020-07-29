#! /usr/bin/python
# coding: utf-8

"""This module is for handling variable meta information.


Example:
    Create an instance of :class:`RemoVariable`::

        from pyremo2.variable import RemoVariable

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

from .codes import get_dict


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

    def __init__(self, id, **kwargs):
        """Creates an instance.

        Args:
            id (int or str): Variable identifier.
        """
        _Variable.__init__(self, **kwargs)
        self._get_info(id)

    def _get_info(self, id):
        self.varinfo = get_dict(id)
        if self.varinfo:
            for prop, value in self.varinfo.items():
                setattr(self, prop, value)
        else:
            raise Exception("Undefined Identifier: {}!".format(id))


def variable(id):
    """Top level function to create a Remo variable.

    Args:
        id (int or str): Variable identifier.

    """
    return RemoVariable(id)
