# from pyremo.physics import specific_humidity, liquid_water_content

import pyremo.physics as prp
import xarray as xr


zds3 = (3.25 - 0.0)/(3.5 - 0.0)
zds4 = (19.2 - 17.5)/(64.0 - 17.5)
zds5 = (77.55 - 64.0)/(194.5 - 64.0)
zdtfak = 0.00065


# freezing and melting point for
# seaice derivation
frozen = 271.37
melt   = 274.16

# from mo_comphy
T0 = 273.16
R = 2.8705E2
RD = 4.6151E2
WCP = 1.005E3
WLK = 2.501E6
WLF = 0.334E6
WLS = 2.835E6
G = 9.80665
RERd = 6.371229E6
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
RDRd = R/RD
RDDrm1 = RD/R - 1.0
EMRdrd = 1.0 - RDRd


def soil_layers(state):
#    if 'TD' not in state:
    tswem = state['TSW']
    tslem = state['TSL']
    td3ge = state.get('TD3', None)
    td4ge = state.get('TD4', None)
    td5ge = state.get('TD5', None)
    tdge = state.get('TD', None)
    cm.print_datinfo('TSW', tswem)
    cm.print_datinfo('TSL', tslem)
    #cm.print_datinfo('td3ge', td3ge)
    #cm.print_datinfo('td4ge', td4ge)
    #cm.print_datinfo('td5ge', td5ge)
    #cm.print_datinfo('tdge', tdge)
    logging.info('adding soil layers')
    soil = physics.soil_layers(tswem, tslem, sd.bodlib.blaem, sd.bodlib.fibem, sd.bodlib.fibge, td3ge, td4ge, td5ge, tdge )
    for field, data in soil.items():
        cm.print_datinfo(field, data)
    state['TD3'] = soil['td3em']
    state['TD4'] = soil['td4em']
    state['TD5'] = soil['td5em']
    state['TD'] = soil['tdem']
    state['TDCL'] = soil['tdclem']
    return state

    

def water_content(ads):
    p = prp.pressure(ads.PS, ads.akh, ads.bkh)
    qwem = prp.liquid_water_content(ads.T, ads.RF, p)
    qdem = prp.specific_humidity(ads.T, ads.RF, p)
    ads['QD'] = qdem
    ads['QW'] = qwem
    ads['QDBL'] = qdem.isel(lev=-1) # last layer will be used for QDBL
    return ads



def seaice(tswem):
    """Simple derivation of seaice fraction from SST.

    A simple linear approximation is done if the water temperature is
    between the freezing and the melting point.
    """
    freezing_range = melt - frozen
    seaem = xr.zeros_like(tswem)
    seaem = xr.where( tswem > melt, 0.0, seaem)
    seaem = xr.where( tswem < frozen, 1.0, seaem)
    seaem = xr.where( (tswem >= frozen) & (tswem <= melt), (melt - tswem) / freezing_range, seaem )

    return seaem



def tsw(tswge):
    tswem = np.ma.where( tswge < frozen, frozen, tswge)
    return tswem

def tsi(tswge):
    tsiem = np.ma.where( tswge > frozen, frozen, tswge)
    return tsiem
