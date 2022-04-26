"""Compute derived variables for cmorization.
"""


import inspect
from warnings import warn

from .. import physics

# units are incompatible: kg m-2 s-1 and mm

# pr = 142 + 143
#    = APRL + APRC
#    = large scale precipitation + convective precipitation


class derivator:
    """derivator static class.

    The derivator class provides access to the derivation of
    variables for cmorization. For each derived variable, a
    static class method should be provided of the same name as
    the cf variable. The parameters of the method should have the
    REMO variable names that are required to derive the cf variable.
    The derive method can then determine which variables to take
    from the dataset and how to derive the variable automatically.
    """

    @classmethod
    def get_function(cls, cf_varname):
        func = getattr(cls, cf_varname)
        return func

    @classmethod
    def get_params(cls, cf_varname):
        func = getattr(cls, cf_varname)
        return list(inspect.signature(func).parameters)

    @classmethod
    def derive(cls, ds, cf_varname):
        func = cls.get_function(cf_varname)
        params = cls.get_params(cf_varname)
        args = (ds[param] for param in params if param in ds)
        warn("computing {} from {}".format(cf_varname, [p for p in params]))
        return func(*args)

    @staticmethod
    def pr(APRL, APRC):
        res = physics.precipitation_flux(APRL, APRC)
        res.name = "pr"
        res.attrs["units"] = "mm"
        return res

    @staticmethod
    def huss(DEW2, PS):
        res = physics.specific_humidity_from_dewpoint(DEW2, PS)
        res.name = "huss"
        res.attrs["units"] = "kg/kg"
        return res

    @staticmethod
    def hurs(DEW2, TEMP2):
        res = physics.relative_humidity_from_dewpoint(DEW2, TEMP2)
        res.name = "hurs"
        res.attrs["units"] = "%"
        return res

    @staticmethod
    def mrros(RUNOFF, DRAIN):
        res = physics.surface_runoff_flux(RUNOFF, DRAIN)
        res.name = "mrros"
        res.attrs["units"] = "mm"
        return res

    @staticmethod
    def rsds(SRADS, SRADSU):
        res = physics.surface_downwelling_shortwave_flux_in_air(SRADS, SRADSU)
        res.name = "rsds"
        res.attrs["units"] = "W m-2"
        return res

    @staticmethod
    def rlds(TRADS, TRADSU):
        res = physics.surface_downwelling_longwave_flux_in_air(TRADS, TRADSU)
        res.name = "rlds"
        res.attrs["units"] = "W m-2"
        return res

    @staticmethod
    def rsdt(SRAD0, SRAD0U):
        res = physics.toa_incoming_shortwave_flux(SRAD0, SRAD0U)
        res.name = "rsdt"
        res.attrs["units"] = "W m-2"
        return res

    @staticmethod
    def evspsbl(EVAP):
        res = physics.water_evapotranspiration_flux(EVAP)
        res.name = "evspsbl"
        res.attrs["units"] = "mm"
        return res
