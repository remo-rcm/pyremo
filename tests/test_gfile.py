import pytest

from pyremo.preproc import cf
from pyremo.tutorial import load_dataset, mpi_esm_tos


@pytest.fixture
def mpi_esm_sources():
    files = {
        "ta": "ta_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "ua": "ua_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "va": "va_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "hus": "hus_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "ps": "ps_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197901010600",
        "orog": "orog_fx_MPI-ESM1-2-HR_historical_r1i1p1f1_gn",
        "sftlf": "sftlf_fx_MPI-ESM1-2-HR_historical_r1i1p1f1_gn",
    }

    file_dict = {k: load_dataset(f).encoding["source"] for k, f in files.items()}
    file_dict["tos"] = mpi_esm_tos().encoding["source"]
    return file_dict


def test_gfile_mpi(mpi_esm_sources):
    gfile = cf.get_gfile(**mpi_esm_sources)
    gds = gfile.gfile("1979-01-01T06:00:00")
    vars = mpi_esm_sources.keys()
    for var in vars:
        assert var in gds
