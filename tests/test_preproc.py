
import pytest
import pyremo as pr
from pyremo.preproc.double_nesting import process_file

from . import requires_pyintorg

@pytest.fixture
def ds():
    return pr.tutorial.load_dataset()


@requires_pyintorg
def test_double_nesting(ds):
    em = pr.domain_info('EUR-44')
    hm = pr.domain_info('EUR-11')
    hm_ds = pr.remo_domain('EUR-11')
    em_ds = pr.remo_domain('EUR-44')
    vc = pr.vc.tables['vc_27lev']
    surflib = pr.data.surflib('EUR-11', crop=False)
    ads = process_file(ds, em, hm, vc, surflib)
    print(ads)
