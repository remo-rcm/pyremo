# -*- coding: utf-8 -*-
# flake8: noqa
import pytest

import pyremo as pr


def test_domain_dataset():

    eur11 = pr.remo_domain("EUR-11")
    eur22 = pr.remo_domain("EUR-22")
    eur44 = pr.remo_domain("EUR-44")
