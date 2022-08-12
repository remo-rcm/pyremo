import os

import pytest
import xarray as xr

import pyremo as pr
from pyremo.core import cli as add_vars_cli
from pyremo.prsint import cli as prsint_cli

from . import requires_pydruint


@pytest.fixture
def temp2():
    ds = pr.tutorial.open_dataset("remo_EUR-11_TEMP2_1hr")
    return ds.encoding["source"]


@pytest.fixture
def tfile():
    ds = pr.tutorial.open_dataset()
    return ds.encoding["source"]


@requires_pydruint
def test_entrypoint_prsint():
    exit_status = os.system("prsint --help")
    assert exit_status == 0


def test_entrypoint_add_vars():
    exit_status = os.system("pradd-vars --help")
    assert exit_status == 0


@requires_pydruint
@pytest.mark.parametrize(
    "split, expected",
    [
        (
            True,
            [
                "e056000p_c130_0100_200601.nc",
                "e056000p_c130_0100_200601.nc",
                "e056000p_c130_0100_200601.nc",
            ],
        ),
        (False, ["e056000p_c130_200601.nc"]),
    ],
)
def test_prsint_cli(tfile, split, expected):
    parser = prsint_cli.prsint_parser()
    args = (tfile, "-id", "056000", "-v", "T", "--plev", "100", "200", "300")
    if split is True:
        args = args + ("-s",)
    args = parser.parse_args(args)
    prsint_cli.prsint(args)
    for f in expected:
        assert os.path.exists(f)


def test_add_vars_cli(temp2):
    parser = add_vars_cli.replace_parser()
    pr.data.surflib("EUR-11").to_netcdf("surflib.nc")
    print(temp2)
    args = (temp2, "surflib.nc", "-v", "BLA", "FIB")
    args = parser.parse_args(args)
    fname = add_vars_cli.replace_variables(args)
    ds = xr.open_dataset(fname)
    assert "BLA" in ds
    assert "FIB" in ds
