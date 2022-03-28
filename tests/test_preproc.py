
import pytest
import pyremo as pr

from . import requires_pyintorg

@pytest.fixture
def ds():
    return pr.tutorial.load_dataset()


@requires_pyintorg
def test_double_nesting(ds):
    print(ds)
