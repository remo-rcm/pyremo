# from pyremo.physics import specific_humidity, liquid_water_content

import numpy as np
import xarray as xr

import pyremo.physics as prp

from . import constants as const

# Constants
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


# def calculate_dp(ps, ak, bk):
#     return ps - (ak[-1] + bk[-1] * ps)


def calculate_zdts(fibem, fibge):
    """height correction"""
    return (fibem - fibge) * zdtfak


def adapt_soil_temperatures(
    tdge, tswge, tslge, td3ge, td4ge, td5ge, fibem, fibge, blaem
):
    # return tdge, tswge, tslge, td3ge, td4ge, td5ge, fibem, fibge, blaem
    iland = (blaem * 10000).astype(int)
    zdts = calculate_zdts(fibem, fibge)
    tslem = xr.where(iland == 0, tswge, tslge - zdts)
    tslem.name = "TSL"
    tswem = tswge.copy()
    tswem.name = "TSW"
    tsnem = xr.where(iland == 0, tswge, tslem)
    tsnem.name = "TSN"
    td3em = xr.where(iland == 0, tswge, tslge - ((tslge - td3ge) * zds3) - zdts)
    td3em.name = "TD3"
    td4em = xr.where(iland == 0, tswge, td4ge - ((td4ge - td5ge) * zds4) - zdts)
    td4em.name = "TD4"
    td5em = xr.where(iland == 0, tswge, td5ge - ((td5ge - tdge) * zds5) - zdts)
    td5em.name = "TD5"
    tdem = td5em.copy()
    tdem.name = "TD"
    tdclem = td5em.copy()
    tdclem.name = "TDCL"
    return tslem, tswem, tsnem, td3em, td4em, td5em, tdem, tdclem


def calculate_seaice(tswem):
    seaem = xr.zeros_like(tswem)
    seaem = xr.where(tswem > melt, 0.0, seaem)
    seaem = xr.where(tswem < frozen, 1.0, seaem)
    seaem = xr.where(
        (tswem >= frozen) & (tswem <= melt), (melt - tswem) / (melt - frozen), seaem
    )
    return seaem


def calculate_tsi(tswem):
    tsiem = tswem.copy()
    tswem = xr.where(tswem < frozen, frozen, tswem)
    tsiem = xr.where(tsiem > frozen, frozen, tsiem)
    return tswem, tsiem


def diagnose_ice_mask(snem):
    return xr.where(snem > 9.50, 1.0, 0.0)


def bodfld(
    psge,
    psem,
    akgm,
    bkgm,
    akem,
    bkem,
    snge,
    wsge,
    wsmxem,
    wlge,
    blaem,
    fibem,
    fibge,
    tswge,
    tslge,
    td3ge,
    td4ge,
    td5ge,
    tdge,
    prgpk=False,
    ipr=None,
    jpr=None,
):
    # Calculate dpge and dpem
    # dpge = calculate_dp(psge, akgm, bkgm)
    # dpem = calculate_dp(psem, akem, bkem)

    # fak1 = 0.07_DP
    # fak2 = 0.21_DP
    # fak3 = 0.72_DP

    # WSGm(ij) = (fak1*zwb1(ij)) + (fak2*zwb2(ij)) + (fak3*zwb3(ij))

    # Initialize arrays
    snem = snge.copy()
    wsem = xr.ufuncs.minimum(wsmxem, wsge * wsmxem)
    wlem = wlge.copy()
    iland = xr.ufuncs.round(blaem * 1.0).astype(
        int
    )  # Assuming FAKINF is 1.0 for simplicity

    # Initialize temperatures
    tslem, tswem, tsnem, td3em, td4em, td5em, tdem, tdclem = adapt_soil_temperatures(
        tswge, tslge, td3ge, td4ge, td5ge, fibem, fibge, iland
    )

    # Calculate sea ice and tsi
    seaem = calculate_seaice(tswem)
    tswem, tsiem = calculate_tsi(tswem)

    # Diagnose ice mask from snow height
    glacem = diagnose_ice_mask(snem)

    return (
        tslem,
        tswem,
        tsnem,
        td3em,
        td4em,
        td5em,
        tdem,
        tdclem,
        seaem,
        tsiem,
        glacem,
        wsem,
        wlem,
        snem,
    )

    # results = bodfld(psge, psem, akgm, bkgm, akem, bkem, snge, wsge, wsmxem, wlge, blaem, fibem, fibge, tswge, tslge, td3ge, td4ge, td5ge, tdge)
    # print(results)


#  & 'FIB     ',   129    ,    1.0   ,      0.0  ,                              &
#   & 'BLA     ',   172    ,    1.0   ,      0.0  ,                              &
#   & 'PHI     ',   253    ,   57.296 ,      0.0  ,                              &
#   & 'LAM     ',   254    ,   57.296 ,      0.0  ,                              &
#   & 'PS      ',   152    ,    1.0   ,      0.0  ,                              &
#   & 'T       ',   130    ,    1.0   ,      0.0  ,                              &
#   & 'QD      ',   133    ,    1.0   ,      0.0  ,                              &
#   & 'QW      ',   246    ,    1.0   ,      0.0  ,                              &
#   & 'U       ',   131    ,    1.0   ,      0.0  ,                              &
#   & 'V       ',   132    ,    1.0   ,      0.0  /

# DATA(GMNAME(I), GMGBNR(I), GMGBFK(I), GMGBBS(I) , I = 11,20) /                 &
#   & 'WB1     ',    39    ,    1.0   ,      0.0  ,                              &
#   & 'WB2     ',    40    ,    1.0   ,      0.0  ,                              &
#   & 'WB3     ',    41    ,    1.0   ,      0.0  ,                              &
#   & 'SN      ',   141    ,    1.0   ,      0.0  ,                              &
#   & 'WL      ',   198    ,    1.0   ,      0.0  ,                              &
#   & 'TS      ',   235    ,    1.0   ,      0.0  ,                              &
#   & 'TD3     ',   139    ,    1.0   ,      0.0  ,                              &
#   & 'TD4     ',   170    ,    1.0   ,      0.0  ,                              &
#   & 'TD5     ',   183    ,    1.0   ,      0.0  ,                              &
#   & 'TD      ',   236    ,    1.0   ,      0.0  /

# DATA(GMNAME(I), GMGBNR(I), GMGBFK(I), GMGBBS(I) , I = 21,NFGM) /               &
#   & 'TSN     ',   238    ,    1.0   ,      0.0  ,                              &
#   & 'SEAICE  ',    31    ,    1.0   ,      0.0  ,                              &
#   & 'SST     ',    34    ,    1.0   ,      0.0  /


#     CALL HIMBLA(TSGm, TSLge, 'TSL ')
# CALL HIMBLA(TSGm, TSWge, 'TSW ')
# CALL HIMBLA(TD3gm, TD3ge, 'TD3 ')
# CALL HIMBLA(TD4gm, TD4ge, 'TD4 ')
# CALL HIMBLA(TD5gm, TD5ge, 'TD5 ')
# CALL HIMBLA(TDGm, TDGe, 'TD  ')
# CALL HIMBLA(TSNgm, TSNge, 'TSN ')
# !
# !     SKIN-RESERVOIR OF PLANTS HORIZONTAL INTERPOLIEREN
# CALL HIMBLA(WLGm, WLGe, 'WL ')
# !
# !     ZONALEN WIND HORIZONTAL INTERPOLIEREN UND
# !
# !     SCHNEEHOEHE HORIZONTAL INTERPOLIEREN
# CALL HIMBLA(SNGm, SNGe, 'SN ')
# !
# !     WASSERGEHALT (OBERE SCHICHT) HORIZONTAL INTERPOLIEREN
# CALL HIMBLA(WSGm, WSGe, 'WS ')
# !
# !     SEE-EIS-BEDECKUNG HORIZONTAL INTERPOLIEREN
# CALL HIMBLA(SEAgm, SEAem, 'SEAICE')
# !
# !     TEMPERATUR - DIFFERENZ TP - TB HORIZONTAL INTERPOLIEREN
# CALL HIMBLA(DTPbgm, DTPbge, 'DTPB  ')
# !
# !     GEOPOTENTIAL CHECK-FLAECHE HORIZONTAL INTERPOLIEREN
# CALL HIOBLA(FICgm, FICge, 'FIC   ', 1)
# !
# !     DRUCKTENDENZ AM BODEN HORIZONTAL INTERPOLIEREN
# CALL HIOBLA(DPDtgm, DPDtge, 'DPDT  ', 1)
