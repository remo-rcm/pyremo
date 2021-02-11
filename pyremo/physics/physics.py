"""physics module that contains some common Remo physics functions.
"""

import numpy as np

try:
    import xarray as xr

    hasXarray = True
except:
    hasXarray = False

from . import constants as C


# statement function from druint
def E(tx):
    """statement function"""
    return C.B1 * np.exp(C.B2W(tx - C.B3) / (tx - C.B4W))


# statement function from druint
def GQD(ex, px):
    """statement function"""
    return C.R / C.RD * ex / (px - (1.0 - C.R / C.RD) * ex)


# from preprocessor RemapToRemo
#!     STATEMENTFUNKTION FUER SAETTIGUNGSDAMPFDRUCK
# FGEW(tx) = B1*EXP(B2W*(tx - B3)/(tx - B4W))
# FGEE(tx) = B1*EXP(B2E*(tx - B3)/(tx - B4E))
# FGQD(ge, p) = RDRd*ge/(p - EMRdrd*ge)


def fg(tx, b2, b4):
    """statement function dummy.

    **Arguments:**
        *tx:*
            temperature
        *b2:*
            dampfdruck
        *b4:*
            dampfdruck

    **Returns:**
        fg: statement function
    """
    return C.B1 * np.exp(b2 * (tx - C.B3) / (tx - b4))


def fgew(tx):
    """statement function

    **Arguments:**
        *tx:*
            temperature

    **Returns:**
        fgew: statement function
    """
    return fg(tx, C.B2W, C.B4W)


def fgee(tx):
    """statement function

    **Arguments:**
        *tx:*
            temperature

    **Returns:**
        fgee: statement function
    """
    return fg(tx, C.B2E, C.B4E)


def fgqd(ge, p):
    """statement function

    Arguments
    ---------
    ge : input data, array-like
        input data for the function
    p : pressure data, array-like
        the pressure field

    Returns
    -------
    fgqd: statement function result
    """
    return C.RDRd * ge / (p - C.EMRdrd * ge)


def relative_humidity(t, qd, ps, ak, bk):
    """computes relative humidity from temperature, pressure and specific humidty (qd)

    This might be similar to https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_specific_humidity.html

    algorithm from druint.addfeld

    **Arguments:**
        *ps:*
            surface pressure field ([Pa], 2d)
        *t:*
            temperature fields ([K], 3d)
        *qd:*
            specific humidity fields ([kg/kg], 3d)
        *ak, bk:*
            hybrid coordinates (1d)

    **Returns:**
        *relhum:*
            relative humidity ([%],3d)
    """
    p = pressure_at_model_levels(ps, ak, bk)
    return qd / fgqd(fgew(t), p)


def vc_full_level(ak, bk):
    """computes hybrid vertical coordinates at full levels (akh, bkh).

    Uses a simple linear interpolation.

    **Arguments:**
        *ak, bk:*
            arrays of hybrid coordinates at level interfaces.

    """
    return 0.5 * (ak[:-1] + ak[1:]), 0.5 * (bk[:-1] + bk[1:])


def pressure_at_model_level(ps, ak, bk):
    """computes pressure at a certain model level.

    Uses surface pressure and vertical hybrid coordinates.

    **Arguments:**
        *ps:*
            surface pressure
        *ak, bk:*
            hybrid coordinates at a certain level
    """
    return ak + bk * ps


def pressure_at_model_levels(ps, ak, bk):
    """computes pressure at all model levels.

    Uses surface pressure and vertical hybrid coordinates.

    **Arguments:**
        *ps:*
            surface pressure field
        *ak, bk:*
            arrays of hybrid coordinates at level interfaces.
    """
    p = np.zeros(ps.shape + ak.shape, dtype=ps.dtype)
    for k in range(0, p.shape[2]):
        p[:, :, k] = pressure_at_model_level(ps, ak[k], bk[k])
    return p


def pressure_at_model_levels_xarray(ps, ak, bk):
    """computes pressure at all model levels.

    Uses surface pressure and vertical hybrid coordinates.
    """
    p = xr.DataArray(
        dims=("lev",) + ps.dims,
        coords={**{"lev": range(1, ak.shape[0] + 1)}, **ps.coords},
    )
    for k, p_lev in enumerate(p):
        p_lev[:] = pressure_at_model_level(ps, ak[k], bk[k])
    return p


# from RemapToRemo addem
#!     BERECHNUNG DER SPEZIFISCHEN FEUCHTE UND WOLKENWASSER
# DO k = 1, KEEM
#  DO ij = 1, IJ2EM
#    phem = 0.5_DP*(AKEm(k) + AKEm(k + 1) + (BKEm(k) + BKEm(k + 1))*PSEm(ij))
#    IF ( TEM(ij, k)>=B3 ) THEN
#      zgqd = FGQD(FGEW(TEM(ij, k)), phem)
#    ELSE
#      zgqd = FGQD(FGEE(TEM(ij, k)), phem)
#    END IF
#    zqdwem = ARFem(ij, k)*zgqd
#    IF ( ARFem(ij, k)<1.0_DP ) THEN
#      QDEm(ij, k) = zqdwem
#      QWEm(ij, k) = 0.
#    ELSE
#      QDEm(ij, k) = zgqd
#      QWEm(ij, k) = zqdwem - zgqd
#    END IF
#  END DO
#!


def specific_humidity(t, relhum, ps, ak, bk):
    """computes specific humidity (qd).

    implements algorithm from RemapToRemo addem.
    """
    p = pressure_at_model_levels(ps, ak, bk)
    fge = np.where(t >= C.B3, fgew(t), fgee(t))
    zgqd = fgqd(fge, p)
    zqdw = relhum * zgqd
    return np.where(relhum < 1.0, zqdw, zgqd)


def liquid_water_content(t, relhum, ps, ak, bk):
    """computes liquid water content (qw).

    implements algorithm from RemapToRemo addem.
    """
    p = pressure_at_model_levels(ps, ak, bk)
    fge = np.where(t >= C.B3, fgew(t), fgee(t))
    zgqd = fgqd(fge, p)
    zqdw = relhum * zgqd
    return np.where(relhum < 1.0, 0.0, zqdw - zgqd)
