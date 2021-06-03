# -*- coding: utf-8 -*-
# flake8: noqa
import pytest
from pyremo import domain as dm



def test_refine():
    # check if all 0.11 domains are consistent with the 0.44 domains
    # for short_name, domain in dm.domains('CORDEX').items():
    #    name = domain.short_name.split('-')[0]+'-44'
    #    print(name)
    #    assert(dm.domain(name) * 1.0 == domain)

    eur11 = dm.remo_domain("EUR-11")
    eur22 = dm.remo_domain("EUR-22")
    eur44 = dm.remo_domain("EUR-44")
    # assert(eur22 == eur44.refine(2))
    # assert(eur11 == eur44.refine(4))
    ## test simple domain math.
    # assert(eur22 == 0.5 * eur44 )
    # assert(eur11 == 0.25 * eur44 )
    # assert(eur11 == 0.5 * eur22 )
    # assert(eur44 == 4 * eur11 )


def test_write():
    domain = dm.remo_domain("EUR-11")
    domain.to_netcdf("EUR-11.nc")


if __name__ == "__main__":
    test_write()
