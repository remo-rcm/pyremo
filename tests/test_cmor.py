import pytest
import pyremo as pr
import cftime as cfdt
import datetime as dt


def test_cftime():
    assert pr.cmor.to_cftime(dt.datetime(2000, 1, 1, 1)) == cfdt.datetime(2000, 1, 1, 1)
