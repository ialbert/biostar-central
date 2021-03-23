import re
import bleach
import logging
import time
import uuid
from functools import wraps
from itertools import islice, count
from datetime import datetime
from calendar import timegm
from django.utils.timezone import utc


logger = logging.getLogger('engine')


def fixcase(name):
    return name.upper() if len(name) == 1 else name.lower()


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def split(text):
    text = ''.join(text.split())
    rows = text.split('+')
    return rows


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def strip_tags(text):
    "Strip html tags from text"
    text = bleach.clean(text, tags=[], attributes={}, styles=[], strip=True)
    return text


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


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)

