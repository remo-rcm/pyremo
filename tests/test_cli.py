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
    "output, expected",
    [
        (
            "plev",
            [
                "e056000p_c130_0100_200601.nc",
                "e056000p_c130_0100_200601.nc",
                "e056000p_c130_0100_200601.nc",
            ],
        ),
        ("input", ["e056000p_c130_200601.nc"]),
    ],
)
def test_prsint_cli(tfile, output, expected):
    parser = cli.prsint_parser()
    args = parser.parse_args(
        (tfile, "-id", "056000", "-v", "T", "--plev", "100", "200", "300", "-o", output)
    )
    cli.prsint(args)
    for f in expected:
        assert os.path.exists(f)
