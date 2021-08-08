



#from pyremo.physics import specific_humidity, liquid_water_content

import pyremo.physics as prp
import xarray as xr


def liquid_water_content(arf, t, ps, ak, bk):
    return None


def specific_humidity(arf, t, ps, ak, bk):
    return None


def pressure(ak, bk, ps):
    """Returns pressure levels

    Computes pressure levels from hybrid sigma coeffiecients
    using the surface pressure ps.

    """
    return ak + bk * ps


def water_content(ads, ak, bk):
    arfem  = ads.RF.values
    tem = ads.T.values
    psem = ads.PS.values
    #qdem, qwem = physics.water_content(arfem, tem, psem, ak, bk)
    #ak = 0.5 * (ak[:-1]+ak[1:])
    #bk = 0.5 * (bk[:-1]+bk[1:])
    qwem = prp.liquid_water_content(tem, arfem, psem, ak, bk)
    qdem = prp.specific_humidity(tem, arfem, psem, ak, bk)
    #state['QD'] = qdem
    #state['QW'] = qwem
    #state['QDBL'] = qdem[:,:,-1] # last layer will be used for QDBL
    return qdem