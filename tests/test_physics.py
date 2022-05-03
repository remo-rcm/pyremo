import pytest

import pyremo as pr
from pyremo import physics


@pytest.fixture
def ds():
    return pr.tutorial.load_dataset()


def test_pressure(ds):
    p = physics.pressure(ds.PS, ds.hyai, ds.hybi, ptop=ds.hyai[0])
    assert (p.isel(nhyi=0) == ds.hyai[0]).all()
    assert (p.isel(nhyi=-1) == ds.PS).all()
    p = physics.pressure(ds.PS, ds.hyam, ds.hybm)
    assert (p.isel(nhym=0) == ds.hyam[0]).all()
    # assert(p[ == ds.hyai[0])
