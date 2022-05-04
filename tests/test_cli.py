import os

import pytest

import pyremo as pr
from pyremo.prsint import cli

from . import requires_pydruint


@pytest.fixture
def tfile():
    ds = pr.tutorial.open_dataset()
    return ds.encoding["source"]


@requires_pydruint
def test_entrypoint():
    exit_status = os.system("prsint --help")
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
    parser = cli.prsint_parser()
    args = (tfile, "-id", "056000", "-v", "T", "--plev", "100", "200", "300")
    if split is True:
        args = args + ("-s",)
    args = parser.parse_args(args)
    cli.prsint(args)
    for f in expected:
        assert os.path.exists(f)
