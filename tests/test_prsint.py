import numpy as np
import pytest

import pyremo as pr

from . import requires_pydruint


@pytest.fixture
def tfile():
    return pr.tutorial.load_dataset()


@requires_pydruint
def test_prsint(tfile):
    from pyremo.prsint import pressure_interpolation

    # define pressure levels in hPa
    plev = [100, 200, 500, 850, 950]
    t_plev = pressure_interpolation(
        tfile.T,
        plev=plev,
        t=tfile.T,
        ps=tfile.PS,
        orog=tfile.FIB,
        a=tfile.hyai,
        b=tfile.hybi,
        keep_attrs=True,
    )
    np.testing.assert_array_equal(np.array(plev) * 100.0, t_plev.plev)
