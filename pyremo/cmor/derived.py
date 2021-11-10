"""Compute derived variables for cmorization.
"""


import numpy as np

# units are incompatible: kg m-2 s-1 and mm

# pr = 142 + 143
#    = APRL + APRC
#    = large scale precipitation + convective precipitation


def mm_to_kg(da):
    """1 kg/m2/s = 86400 mm/day."""
    return da * 1.0 / 86400.0


def pr(aprl, aprc):
    """Precipitation flux"""
    return aprl + aprc


# remo: dew = dew2 (168)
def water_vapour(t):
    """Partial pressure of water vapour e [Pa].

    Computes water vapour from  temperature.

    references: https://doi.org/10.1175/1520-0450(1967)006%3C0203:OTCOSV%3E2.0.CO;2

    compare to: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_vapor_pressure.html#metpy.calc.saturation_vapor_pressure

    """
    T_0 = 273.15
    T_rw = 35.86  # over water
    a = 17.269
    # cdo -mulc,610.78 -exp -div -mulc,17.5 -subc,273.15 a
    return 610.78 * np.exp(a * (t - T_0) / (t - T_rw))


def specific_humidity(dew, ps):
    """Specific humidity `huss` [kg/kg].

    Computes specific humidity `huss` from dewpoint temperature and pressure.

    reference: https://de.wikipedia.org/wiki/Luftfeuchtigkeit#Spezifische_Luftfeuchtigkeit

    compare to: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.specific_humidity_from_dewpoint.html

    """
    e = water_vapour(dew)
    return (0.622 * e) / (ps - 0.378 * e)


def relative_humidity(dew, t2m):
    """Relative humidity `hurs` [%].

    Computes relative humidity `hurs` from dewpoint and air temperature.

    compare to: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_dewpoint.html

    """
    e_dew = water_vapour(dew)
    e_t2m = water_vapour(t2m)
    return e_dew / e_t2m


# old cdo formulas
# from helper
#
# if var == 'pr':
#    f1 = os.path.join(store_dir,file_template % 142); f2 = os.path.join(store_dir, file_template % 143)
#    lg.info(shell('cdo add %s %s %s' % (f1, f2, file_out)))
# elif var == 'huss':
#    # compute huss with formula from philip's script
#    f1 = os.path.join(store_dir,file_template % 168); f2 = os.path.join(store_dir, file_template % 134)
#    lg.info(shell('cdo -mulc,610.78 -exp -div -mulc,17.5 -subc,273.15 %s -subc,35.86 %s %s' % (f1,
#        f1, os.path.join(tmp_dir,'EDDEW'))))
#    lg.info(shell('cdo div -mulc,0.622 ' + os.path.join(tmp_dir,'EDDEW')+' -sub ' + f2 + ' -mulc,0.378 ' +
#        os.path.join(tmp_dir,'EDDEW') + ' ' + file_out))
# elif var == 'hurs':
#    # compute hurs with formula from philip's script
#    f1 = os.path.join(store_dir, file_template % 168); f2 = os.path.join(store_dir,file_template % 167)
#    lg.info(shell('cdo %s -mulc,610.78 -exp -div -mulc,17.5 -subc,273.15 %s -subc,35.86 %s %s' % (omp_flag, f1,
#        f1, os.path.join(tmp_dir,'EDDEW'))))
#    lg.info(shell('cdo %s -mulc,610.78 -exp -div -mulc,17.5 -subc,273.15 %s -subc,35.86 %s %s' % (omp_flag, f2, f2, os.path.join(tmp_dir,'EDT2M'))))
#    lg.info(shell('cdo %s div ' % (omp_flag) + os.path.join(tmp_dir,'EDDEW')+' '+ os.path.join(tmp_dir,'EDT2M')+' ' + file_out))
## TODO: here we can use the codes from the pre list...
# elif var == 'mrros':
#    f1 = os.path.join(store_dir,file_template % 160); f2 = os.path.join(store_dir,file_template % 53)
#    lg.info(shell('cdo sub %s %s %s' % (f1, f2, file_out)))
# elif var == 'rsds':
#    f1 = os.path.join(store_dir, file_template % 176); f2 = os.path.join(store_dir, file_template % 204)
#    lg.info(shell('cdo sub %s %s %s' % (f1, f2, file_out)))
# elif var == 'rlds':
#    f1 = os.path.join(store_dir, file_template % 177); f2 = os.path.join(store_dir, file_template % 205)
#    lg.info(shell('cdo sub %s %s %s' % (f1, f2, file_out)))
# elif var == 'rsdt':
#    f1 = os.path.join(store_dir, file_template % 178); f2 = os.path.join(store_dir, file_template % 203)
#    lg.info(shell('cdo sub %s %s %s' % (f1, f2, file_out)))
# elif var == 'evspsbl':
#    f1 = os.path.join(store_dir,  file_template % 182)
#    lg.info(shell('cdo mulc,-1 %s %s' % (f1, file_out)))
