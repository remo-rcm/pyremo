"""this module handles absolute calendars that occur in REMO datasets.
"""

import datetime as dt
import math


# roundTime from here:
# https://stackoverflow.com/questions/3463930/how-to-round-the-minute-of-a-datetime-object-python
def roundTime(datetime=None, roundTo=60):
    """Round a datetime object to any time lapse in seconds
    dt : datetime.datetime object, default now.
    roundTo : Closest number of seconds to round to, default 1 minute.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me.
    """
    if datetime is None:
        datetime = dt.datetime.now()
    seconds = (datetime.replace(tzinfo=None) - datetime.min).seconds
    rounding = (seconds + roundTo / 2) // roundTo * roundTo
    return datetime + dt.timedelta(0, rounding - seconds, -datetime.microsecond)


class AbsoluteCalendar:
    """Absolute calendar to handle absolute dates.

    Forcing Files will probably contain absolute dates. This is not cf-standard so
    we handle this manually (not using the netCDF4 num2date/date2num, those work only
    for relative calendars.).
    """

    def __init__(self, units=None, calendar=None):
        if units is None:
            self.units = "day as %Y%m%d.%f"
        if calendar is None:
            self.calendar = "proleptic_gregorian"
        self.fmt = "%Y%m%d"  # self.units.split()[2]

    def ncattrs_dict(self):
        return {
            "standard_name": "time",
            "units": self.units,
            "calendar": self.calendar,
            "axis": "T",
        }

    def date2num(self, datetime):
        """convert a datetime object to an absolute numeric date value."""
        delta = dt.timedelta(
            hours=datetime.hour, minutes=datetime.minute, seconds=datetime.second
        )
        frac = delta.total_seconds() / dt.timedelta(days=1).total_seconds()
        return float(datetime.strftime(self.fmt)) + frac

    def num2date(self, num, use_cftime=False, roundTo=60, calendar="standard"):
        """convert a numeric absolute date value to a datetime object."""
        frac, whole = math.modf(num)
        date_str = str(int(whole))
        if date_str[6:8] == "00":
            date_str = date_str[0:6] + "15"
        # date = pd.to_datetime(date_str, format=self.fmt) #dt.datetime.strptime(date_str, self.fmt)
        date = dt.datetime.strptime(date_str[0:8], self.fmt)
        datetime = roundTime(
            date + dt.timedelta(seconds=dt.timedelta(days=1).total_seconds() * frac),
            roundTo=roundTo,
        )
        if use_cftime:
            import cftime

            datetime = cftime.datetime(
                datetime.year,
                datetime.month,
                datetime.day,
                datetime.hour,
                datetime.minute,
                calendar=calendar,
            )
        return datetime


def parse_dates(ds, use_cftime=False, calendar="standard"):
    """Update the time axis of a REMO dataset.

    Updates the time axis of a REMO dataset containing an absolute time axis.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with absolute time axis.
    use_cftime: bool
        Use cftime objects instead of datetime objects.
    calendar: str
        CF calendar if use_cftime is True.
    """
    ds = ds.copy()
    ds["time"] = parse_absolute_time(ds.time, use_cftime=use_cftime, calendar=calendar)
    return ds


def parse_absolute_time(time, use_cftime=False, calendar="standard"):
    """Update a time axis containg fractional absolute dates.

    Updates fractional absolute dates to relative dates.

    Parameters
    ----------
    time : array like or xr.DataArray
        Time axis containing absolute dates as float or int.
    use_cftime: bool
        Use cftime objects instead of datetime objects.
    calendar: str
        CF calendar if use_cftime is True.
    """
    parser = AbsoluteCalendar()
    return [parser.num2date(date, use_cftime, calendar=calendar) for date in time]
