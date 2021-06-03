"""this module handles absolute calendars (for a-files.)
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
   if datetime == None : datetime = dt.datetime.now()
   seconds = (datetime.replace(tzinfo=None) - datetime.min).seconds
   rounding = (seconds+roundTo/2) // roundTo * roundTo
   return datetime + dt.timedelta(0,rounding-seconds,-datetime.microsecond)



class AbsoluteCalendar():
    """Absolute calendar to handle absolute dates.

    Forcing Files will probably contain absolute dates. This is not cf-standard so
    we handle this manually (not using the netCDF4 num2date/date2num, those work only
    for relative calendars.).
    """

    def __init__(self, units=None, calendar=None):
        if units is None:
            self.units = 'day as %Y%m%d.%f'
        if calendar is None:
            self.calendar = 'proleptic_gregorian'
        self.fmt = '%Y%m%d'#self.units.split()[2]

    def ncattrs_dict(self):
        return {'standard_name': 'time', 'units': self.units,
                'calendar':self.calendar, 'axis': 'T'}

    def date2num(self, datetime):
        """convert a datetime object to an absolute numeric date value.
        """
        delta = dt.timedelta(hours=datetime.hour, minutes=datetime.minute,
                seconds=datetime.second)
        frac = delta.total_seconds()/dt.timedelta(days=1).total_seconds()
        return float(datetime.strftime(self.fmt)) + frac

    def num2date(self, num):
        """convert a numeric absolute date value to a datetime object.
        """
        frac, whole = math.modf(num)
        date_str = str(int(whole))
        #date = pd.to_datetime(date_str, format=self.fmt) #dt.datetime.strptime(date_str, self.fmt)
        date = dt.datetime.strptime(date_str, self.fmt)
        return roundTime(date + dt.timedelta(seconds = dt.timedelta(days=1).total_seconds() * frac))


