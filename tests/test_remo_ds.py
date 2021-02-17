# -*- coding: utf-8 -*-
from pyremo import remo_ds as rds

test_data = './tests/data/e056111t2006010100'


def test_remo_ds():
    filename = test_data
    ds = rds.open_remo_dataset(filename, update_meta=True)
    assert 'T' in ds
    assert ds.T.code == 130
