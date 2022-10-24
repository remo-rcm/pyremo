# -*- coding: utf-8 -*-

import pyremo as pr


def test_remo_ds():
    # pytest.importorskip('cdo')
    ds = pr.tutorial.open_dataset()
    assert "T" in ds
    assert ds.T.code == 130


def test_preprocess():
    ds = pr.remo_domain("EUR-44", dummy=True).rename(dummy="T")
    ds = pr.preprocess(ds)
    attrs = pr.codes.get_dict("T")
    for key, value in attrs.items():
        if value is not None:
            assert ds.T.attrs[key] == attrs[key]
