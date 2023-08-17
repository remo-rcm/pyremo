import pint

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

# register special remo units
ureg.define("percent = 0.01*count = %")


def units(units):
    if units == "%":
        return ureg.percent
    elif isinstance(units, str) and hasattr(ureg, units):
        return getattr(ureg, units)
    else:
        print("Warning, unknown unit: {}".format(units))
        return ureg.dimensionless


cf_units = {"kelvin": "K", "celsius": "degC"}


# legs1 = Q_(np.asarray([3., 4.]), 'meter')
# print(legs1)
# print(legs1.dimensionality)
