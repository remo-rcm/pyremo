import pytest
import xarray as xr
import pyremo as pr
from pyremo import cmor as prcmor
import cftime as cfdt
import datetime as dt
import cordex as cx
from cordex.tables import cordex_cmor_table, cmip6_cmor_table


def test_cftime():
    assert prcmor.to_cftime(dt.datetime(2000, 1, 1, 1)) == cfdt.datetime(2000, 1, 1, 1)


def test_cmorizer_fx():
    ds = pr.data.surflib("EUR-11")
    eur11 = cx.cordex_domain("EUR-11")
    ds = ds.assign_coords({"lon": eur11.lon, "lat": eur11.lat})
    filename = prcmor.cmorize_variable(
        ds,
        "orog",
        cmor_table=cmip6_cmor_table("CMIP6_fx"),
        dataset_table=cordex_cmor_table("CORDEX_remo_example"),
        CORDEX_domain="EUR-11",
        time_units=None,
        allow_units_convert=True,
    )
    output = xr.open_dataset(filename)
    assert 'orog' in output
    
    
def test_cmorizer_mon():
    ds = pr.tutorial.open_dataset("remo_EUR-11_TEMP2_mon")
    eur11 = cx.cordex_domain("EUR-11")
    ds = ds.assign_coords({"lon": eur11.lon, "lat": eur11.lat})
    filename = prcmor.cmorize_variable(
        ds,
        "tas",
        cmor_table=cmip6_cmor_table("CMIP6_Amon"),
        dataset_table=cordex_cmor_table("CORDEX_remo_example"),
        CORDEX_domain="EUR-11",
        time_units=None,
        allow_units_convert=True,
    )
    output = xr.open_dataset(filename)
    assert 'tas' in output


@pytest.mark.parametrize('table', ["CMIP6_day", "CMIP6_3hr"])
def test_cmorizer_subdaily(table):
    ds = pr.tutorial.open_dataset("remo_EUR-11_TEMP2_1hr")
    eur11 = cx.cordex_domain("EUR-11")
    ds = ds.assign_coords({"lon": eur11.lon, "lat": eur11.lat})
    filename = prcmor.cmorize_variable(
        ds,
        "tas",
        cmor_table=cmip6_cmor_table(table),
        dataset_table=cordex_cmor_table("CORDEX_remo_example"),
        CORDEX_domain="EUR-11",
        time_units=None,
        allow_units_convert=True,
    )
    output = xr.open_dataset(filename)
    assert 'tas' in output