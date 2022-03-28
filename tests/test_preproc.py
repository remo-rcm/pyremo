
import pytest

import pyremo as pr

@pytest.fixture
def ds():
    return pr.tutorial.load_dataset()


def test_double_nesting(ds):
    print(ds)
