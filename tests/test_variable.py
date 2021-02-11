#!/usr/bin/env python

"""Tests for `pyremo` package."""

import pytest


from pyremo import variable as var


test_data = './data/e056111t2006010100'


def test_variable():
    """Test variable creation"""
    temp = var.RemoVariable('T')
    temp1 = var.RemoVariable(130)
    assert (temp.code == temp1.code)
    assert (temp.variable == temp1.variable)
    temp2 = var.variable('T')
    assert (temp.variable == temp2.variable)



def test_data_variable_ieg():
    """Test a data variable"""
    pytest.importorskip('cdo')
    from cdo import Cdo
    ds = Cdo().copy(options='-t remo', input=test_data, returnXDataset=True)
    temp = var.RemoVariable('T', ds.T)
