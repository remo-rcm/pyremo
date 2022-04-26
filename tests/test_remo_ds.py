# -*- coding: utf-8 -*-

import pyremo as pr


def test_remo_ds():
    # pytest.importorskip('cdo')
    ds = pr.tutorial.open_dataset()
    assert "T" in ds
    assert ds.T.code == 130
