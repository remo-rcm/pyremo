import pytest

import pyremo as pr
from pyremo.preproc import gfile, remap
from pyremo.tutorial import mpi_esm, mpi_esm_tos

from . import requires_pyintorg


@pytest.fixture
def ds():
    return pr.tutorial.load_dataset()


@pytest.fixture
def gcm_data():
    return gfile(mpi_esm(use_cftime=True), tos=mpi_esm_tos(use_cftime=True).tos)


@requires_pyintorg
def test_double_nesting(ds):
    from pyremo.preproc.double_nesting import process_file

    em = pr.domain_info("EUR-44")
    hm = pr.domain_info("EUR-11")
    vc = pr.vc.tables["vc_27lev"]
    surflib = pr.data.surflib("EUR-11", crop=False)
    process_file(ds, em, hm, vc, surflib, write=False)


@requires_pyintorg
def test_cf_preproc(gcm_data):
    domain_info = pr.domain_info("EUR-11")
    surflib = pr.data.surflib("EUR-11", crop=False)
    vc = pr.vc.tables["vc_27lev"]
    ads = remap(gcm_data, domain_info, vc, surflib)
    assert "T" in ads
