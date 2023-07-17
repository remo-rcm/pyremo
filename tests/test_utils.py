import datetime as dt

import pyremo as pr
from pyremo import cal


def test_absolute_calendar():
    calendar = cal.AbsoluteCalendar()
    assert dt.datetime(1979, 1, 1) == calendar.num2date(19790101)
    assert dt.datetime(1979, 1, 1, 12) == calendar.num2date(19790101.5)
    date = dt.datetime(1979, 1, 1, 12)
    assert calendar.date2num(date) == 19790101.5
    date = dt.datetime(1979, 1, 1)
    assert calendar.date2num(date) == 19790101
    assert [
        dt.datetime(1979, 1, 1, 0, 0),
        dt.datetime(1979, 1, 1, 12, 0),
    ] == pr.parse_absolute_time([19790101, 19790101.5])
