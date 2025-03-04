import pytest

import pandas as pd
import pyremo as pr
from pyremo.preproc import get_gcm_dataset, remap
from pyremo.tutorial import load_dataset, mpi_esm, mpi_esm_tos
from pyremo.preproc.utils import get_filename

from . import requires_pyintorg


@pytest.fixture
def ds():
    return pr.tutorial.load_dataset()


# @pytest.fixture
def gcm_ds1():
    return get_gcm_dataset(
        mpi_esm(use_cftime=True), tos=mpi_esm_tos(use_cftime=True).tos
    )


# @pytest.fixture
def gcm_ds2():
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
    return get_gcm_dataset(file_dict, tos=mpi_esm_tos(use_cftime=True).tos)


# @requires_pyintorg
# def test_double_nesting(ds):
#    from pyremo.preproc.double_nesting import process_file
#
#    em = pr.domain_info("EUR-44")
#    hm = pr.domain_info("EUR-11")
#    vc = pr.vc.tables["vc_27lev"]
#    surflib = pr.data.surflib("EUR-11", crop=False)
#    process_file(ds, em, hm, vc, surflib, write=False)


@requires_pyintorg
@pytest.mark.parametrize("gcm_ds", [gcm_ds1(), gcm_ds2()])
def test_cf_preproc(gcm_ds):
    domain_info = pr.domain_info("EUR-11")
    surflib = pr.data.surflib("EUR-11", crop=False)
    vc = pr.vc.tables["vc_27lev"]
    ads = remap(gcm_ds, domain_info, vc, surflib)
    assert "T" in ads


def test_filename():
    assert (
        get_filename(expid="000000", date=pd.Timestamp("1979-01-01"))
        == "a000000a1979010100.nc"
    )
