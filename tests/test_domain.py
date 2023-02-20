# -*- coding: utf-8 -*-
# flake8: noqa
import pytest

import pyremo as pr


def test_domain_dataset():
    eur11 = pr.remo_domain("EUR-11")
    eur22 = pr.remo_domain("EUR-22")
    eur44 = pr.remo_domain("EUR-44")


def test_magic():
    assert pr.magic_numbers(20) == [
        3,
        5,
        7,
        9,
        11,
        13,
        17,
        19,
        21,
        25,
        31,
        33,
        37,
        41,
        49,
        51,
        55,
        61,
        65,
        73,
    ]
