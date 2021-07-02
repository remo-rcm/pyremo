# -*- coding: utf-8 -*-
import pyremo as pr

import pytest

test_data = './tests/data/e056111t2006010100'


def test_remo_ds():
    pytest.importorskip('cdo')
    filename = test_data
    ds = pr.open_remo_dataset(filename, update_meta=True)
    assert 'T' in ds
    assert ds.T.code == 130
