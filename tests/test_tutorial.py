"""Test loading of tutorial data
"""

import pyremo as pr


def test_mpi_esm():
    mpi_esm = pr.tutorial.mpi_esm()
    for var in [
        "ta",
        "ua",
        "va",
        "ps",
        "hus",
        "orog",
        "sftlf",
    ]:
        assert var in mpi_esm


def test_remo_example():
    remo = pr.tutorial.open_dataset("remo_EUR-44")
    for var in ["T", "PS", "U", "V", "FIB", "BLA", "QD"]:
        assert var in remo
