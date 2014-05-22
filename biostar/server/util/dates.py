"""
Utility module for date format conversion.
"""
from calendar import timegm
from datetime import datetime, timedelta
from os import environ

from biostar.apps.posts.models import Post
from .exceptions import NoDayZeroError


def datetime_to_iso(date):
    """
    Converts a datetime to the ISO8601 format, like: 2014-05-20T06:11:41.733900.

    Parameters:
    date -- a `datetime` instance.
    """
    if not isinstance(date, datetime):
        date = datetime.combine(date, datetime.min.time())
    return date.isoformat()


def datetime_to_unix(date):
    """
    Converts a datetime to a Unix timestamp , like: 1400566301.

    Parameters:
    date -- a `datetime` instance.
    """
    return timegm(date.timetuple())


def unix_to_datetime(timestamp):
    """
    Converts a Unix timestamp (like: 1400566301) to a datetime.

    Parameters:
    timestamp -- a Unix timestamp like: 1400566301.
    """
    return datetime.fromtimestamp(float(timestamp))


def days_after_day_zero_to_datetime(days):
    """
    Converts a date expressed as number of days after day-0 (the date of the first ever post) to
    `datetime`.

    Params:
    days -- number of days after day-0 (the date of the first post ever).
    """
    day_zero = _find_day_zero()
    date = day_zero + timedelta(days=int(days))
    return date


def _find_day_zero():
    """
    Finds the day-0, which is the date of the first post ever.
    The day-0 is computer only one time and then cached in an environment var.
    """
    day_zero = environ.get('ZERO_DAY', None)
    if day_zero:
        return unix_to_datetime(day_zero)

    try:
        first_post = Post.objects.order_by('creation_date').only('creation_date')[0]
        day_zero = first_post.creation_date
    except IndexError:
        # In Python 3: raise NoDayZeroError('...') from IndexError
        raise NoDayZeroError('No posts yet.')

    environ['ZERO_DAY'] = '{}'.format(datetime_to_unix(day_zero))
    return day_zero