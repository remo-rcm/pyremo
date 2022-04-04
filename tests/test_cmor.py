import pytest
import pyremo as pr
from pyremo import cmor as prcmor
import cftime as cfdt
import datetime as dt
import cordex as cx


def test_cftime():
    assert prcmor.to_cftime(dt.datetime(2000, 1, 1, 1)) == cfdt.datetime(2000, 1, 1, 1)


def test_cmorizer():
    surflib = pr.data.surflib("EUR-11")
    eur11 = cx.cordex_domain("EUR-11")
    surflib = surflib.assign_coords({"lon": eur11.lon, "lat": eur11.lat})
    filename = prcmor.cmorize_variable(
        surflib,
        "orog",
        "fx",
        cx.cordex_cmor_table("remo_example"),
        CORDEX_domain="EUR-11",
        time_units=None,
        allow_units_convert=True,
    )
