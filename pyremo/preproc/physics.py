# from pyremo.physics import specific_humidity, liquid_water_content

import numpy as np
import xarray as xr

import pyremo.physics as prp

from . import constants as const

zds3 = (3.25 - 0.0) / (3.5 - 0.0)
zds4 = (19.2 - 17.5) / (64.0 - 17.5)
zds5 = (77.55 - 64.0) / (194.5 - 64.0)
zdtfak = 0.00065


# freezing and melting point for
# seaice derivation
frozen = 271.37
melt = 274.16

# from mo_comphy
T0 = 273.16
R = 2.8705e2
RD = 4.6151e2
WCP = 1.005e3
WLK = 2.501e6
WLF = 0.334e6
WLS = 2.835e6
G = 9.80665
RERd = 6.371229e6
STAg = 86164.09054


# from mo_compar
B1 = 610.78
B2W = 17.2693882
B2E = 21.8745584
B3 = 273.16
B4E = 7.66
B4W = 35.86

# from mo_comhig
PI = 3.141592654
RADdeg = 57.29577951
DEGrad = 0.0174532925
RDRd = R / RD
RDDrm1 = RD / R - 1.0
EMRdrd = 1.0 - RDRd


# def soil_layers(state):
# if 'TD' not in state:
#    tswem = state['TSW']
#    tslem = state['TSL']
#    td3ge = state.get('TD3', None)
#    td4ge = state.get('TD4', None)
#    td5ge = state.get('TD5', None)
#    tdge = state.get('TD', None)
#    #cm.print_datinfo('TSW', tswem)
#    #cm.print_datinfo('TSL', tslem)
#    #cm.print_datinfo('td3ge', td3ge)
#    #cm.print_datinfo('td4ge', td4ge)
#    #cm.print_datinfo('td5ge', td5ge)
#    #cm.print_datinfo('tdge', tdge)
#    #logging.info('adding soil layers')
#    ###soil = prp.soil_layers(tswem, tslem, sd.bodlib.blaem, sd.bodlib.fibem, sd.bodlib.fibge, td3ge, td4ge, td5ge, tdge )
#    #for field, data in soil.items():
#    #    cm.print_datinfo(field, data)
#    #state['TD3'] = soil['td3em']
#    #state['TD4'] = soil['td4em']
#    #state['TD5'] = soil['td5em']
#    #state['TD'] = soil['tdem']
#    #state['TDCL'] = soil['tdclem']
#    return state


def water_content(t, rf, ps, ak, bk):
    p = prp.pressure(ps, ak, bk)
    qwem = prp.liquid_water_content(t, rf, p)
    qdem = prp.specific_humidity(t, rf, p)
    qwem.name = "QW"
    qdem.name = "QD"
    # ads['QD'] = qdem
    # ads['QW'] = qwem
    # ads['QDBL'] = qdem.isel(lev=-1) # last layer will be used for QDBL
    qdbl = qdem.isel(lev=-1)
    qdbl.name = "QDBL"
    return xr.merge([qwem, qdem, qdbl])
    # return ads


def seaice(tswem):
    """Simple derivation of seaice fraction from SST.

    A simple linear approximation is done if the water temperature is
    between the freezing and the melting point.
    """
    freezing_range = melt - frozen
    seaem = xr.zeros_like(tswem)
    seaem = xr.where(tswem > melt, 0.0, seaem)
    seaem = xr.where(tswem < frozen, 1.0, seaem)
    seaem = xr.where(
        (tswem >= frozen) & (tswem <= melt), (melt - tswem) / freezing_range, seaem
    )
    seaem.name = "SEAICE"
    return seaem


def tsw(tsw):
    return xr.where(tsw < frozen, frozen, tsw)


def tsi(tsw):
    tsi = xr.where(tsw > frozen, frozen, tsw)
    tsi.name = "TSI"
    return tsi


def _water_content(arfem, tem, psem, akem, bkem):
    """computes qd (specific humidity) and qw (liquid water) from relative humidty (arfem).

    Python implementation of original Fortran source in `addem`.
    """
    # pressure height?!
    # phem = 0.5*(AKEm(k) + AKEm(k + 1) + (BKEm(k) + BKEm(k + 1))*PSEm(ij))
    # IF ( TEM(ij, k)>=B3 ) THEN
    #  zgqd = FGQD(FGEW(TEM(ij, k)), phem)
    # ELSE
    #  zgqd = FGQD(FGEE(TEM(ij, k)), phem)
    # END IF
    # zqdwem = ARFem(ij, k)*zgqd
    # IF ( ARFem(ij, k)<1.0_DP ) THEN
    #  QDEm(ij, k) = zqdwem
    #  QWEm(ij, k) = 0.
    # ELSE
    #  QDEm(ij, k) = zgqd
    #  QWEm(ij, k) = zqdwem - zgqd
    # END IF
    qdem = np.zeros(arfem.shape, dtype=arfem.dtype)
    qwem = np.zeros(arfem.shape, dtype=arfem.dtype)
    # print_data('arfem', arfem)
    # print_data('qdem', qdem)
    # print_data('qwem', qdem)
    for k in range(arfem.shape[2]):
        print(k, B3)
        phem = 0.5 * (akem[k] + akem[k + 1] + (bkem[k] + bkem[k + 1]) * psem)
        # print_data('tem',tem[:,:,k])
        zgqd = np.where(
            tem[:, :, k] >= B3,
            fgqd(fgew(tem[:, :, k]), phem),
            fgqd(fgee(tem[:, :, k]), phem),
        )
        zqdwem = arfem[:, :, k] * zgqd
        # print_data('zqdwem',zqdwem)
        # print_data('arfem',arfem[:,:,k])
        qdem[:, :, k] = np.where(arfem[:, :, k] < 1.0, zqdwem, zgqd)
        qwem[:, :, k] = np.where(arfem[:, :, k] < 1.0, 0.0, zqdwem - zgqd)

    # print_data('qdem', qdem)
    # print_data('qwem', qwem)

    return qdem, qwem


# STATEMENTFUNKTION FUER SAETTIGUNGSDAMPFDRUCK
def fgew(tx):
    """magnus formula"""
    return B1 * np.exp(B2W * (tx - B3) / (tx - B4W))


def fgee(tx):
    """magnus formula"""
    return B1 * np.exp(B2E * (tx - B3) / (tx - B4E))


def fgqd(ge, p):
    """magnus formula"""
    return RDRd * ge / (p - EMRdrd * ge)


def addem_remo(tds):
    # ws auf relative bodenfeucht umrechnen
    wsem = tds.WS.where(tds.WS < 1.0e9, 0.0) / tds.WSMX
    wsem.name = "WS"
    dtpbem = tds.T.isel({const.lev_input: -1}) - tds.TSL
    dtpbem.name = "DTPB"
    return xr.merge([wsem, dtpbem]).squeeze(drop=True)


# def bodfld_remo(ads, surflib):
#    #  RELATIVES WS IN ABSOLUTES WS ZURUECKRECHNEN
#    wshm = ads.WS * surflib.WSMX
#    tslhm = ads.T.isel({const.lev_input: -1}) - dtpbeh(ij) * dphm(ij) / dpeh(ij)


def from_surflib(surflib):
    pass


def derive_soil_temperatures(tds, ads):
    pass


#  dpeh(ij) = pseh(ij) - GETP(akem(KEEM),bkem(KEEM),pseh(ij),akem(1))
#  dphm(ij) = pshm(ij) - GETP(akhm(KEHM),bkhm(KEHM),pshm(ij),akhm(1))


# DO ij = 1 , IJ2HM
#   tswhm(ij) = tsweh(ij)
#   tsihm(ij) = tsieh(ij)
#   tslhm(ij) = thm(ij,KEHM) - dtpbeh(ij)*dphm(ij)/dpeh(ij)
# ENDDO
# !
# DO ij = 1 , IJ2HM
#   zdts(ij) = tslhm(ij) - tsleh(ij)
#   tsnhm(ij) = tsneh(ij) + zdts(ij)
#   td3hm(ij) = td3eh(ij) + zdts(ij)
#   td4hm(ij) = td4eh(ij) + zdts(ij)
#   td5hm(ij) = td5eh(ij) + zdts(ij)
#   tdhm(ij) = tdeh(ij) + zdts(ij)
#   tdclhm(ij) = tdcleh(ij) + zdts(ij)
# ENDDO
