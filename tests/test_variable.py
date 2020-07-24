#!/usr/bin/env python

"""Tests for `pyremo2` package."""

import pytest


from pyremo2 import variable as var



def test_variable():
    """Test variable creation"""
    temp = var.RemoVariable('T')
    temp1 = var.RemoVariable(130)
    assert (temp.code == temp1.code)
    assert (temp.variable == temp1.variable)
