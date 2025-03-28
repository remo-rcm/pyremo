import pandas as pd

from pyremo import codes


def test_codes():
    expected = {
        "variable": "T",
        "description": "temperature",
        "units": "K",
        "positive": None,
        "layer": 110.0,
        "time_cell_method": "point",
        "cf_name": "ta",
        "code": 130,
        "module": "atmos",
    }
    assert codes.get_dict(130) == expected
    assert codes.get_dict("T") == expected


def test_code_search():
    df = pd.DataFrame(
        {
            "code": 167,
            "variable": "TEMP2",
            "description": "2m temperature",
            "units": "K",
            "positive": None,
            "layer": 1.0,
            "time_cell_method": "mean",
            "cf_name": "tas",
            "module": "atmos",
        },
        index=[0],
    )
    assert df.equals(codes.search(code=167).reset_index(drop=True))
