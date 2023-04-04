from pyremo import output_pattern


def test_output_pattern():
    assert "e*e_c*_198001.nc" == output_pattern("e", date="198001")
    assert "e000000e_c130_198001.nc" == output_pattern(
        "e", expid=0, code=130, date="198001"
    )
    assert "e000000e_c*_198001.nc" == output_pattern("e", expid=0, date="198001")
    assert "e*e_c*_1980*.nc" == output_pattern("e", date="1980")
    assert "e065000f196801*.nc" == output_pattern(
        "f", expid=65000, code=1, date="196801"
    )
    assert "e065000g196801*.nc" == output_pattern("g", expid="065000", date="196801")
    assert "e000000t196801*.nc" == output_pattern("t", expid=0, date="196801")
    assert "e000000t1968010100.nc" == output_pattern("t", expid=0, date="1968010100")
    assert "e000000t196801.nc" == output_pattern("t", expid=0, date="196801", end="")
