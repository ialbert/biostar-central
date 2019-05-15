import re
import bleach
import logging
import uuid

from datetime import datetime, timedelta


from django.template import loader
from django.utils.timezone import utc


def fixcase(name):
    return name.upper() if len(name) == 1 else name.lower()


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def split_tags(text):

    capitalize = lambda txt: txt.upper() if len(txt) == 1 else txt
    return [capitalize(x) for x in text.split(",") if len(x)]


def render(name, **kwds):
    "Helper function to render a template"
    tmpl = loader.get_template(name)
    page = tmpl.render(kwds)
    return page


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def strip_tags(text):
    "Strip html tags from text"
    text = bleach.clean(text, tags=[], attributes={}, styles=[], strip=True)
    return text


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)


def time_ago(date):

    if not date:
        return ''
    delta = now() - date
    if delta < timedelta(minutes=1):
        return 'just now'
    elif delta < timedelta(hours=1):
        unit = pluralize(delta.seconds // 60, "minute")
    elif delta < timedelta(days=1):
        unit = pluralize(delta.seconds // 3600, "hour")
    elif delta < timedelta(days=30):
        unit = pluralize(delta.days, "day")
    elif delta < timedelta(days=90):
        unit = pluralize(int(delta.days / 7), "week")
    elif delta < timedelta(days=730):
        unit = pluralize(int(delta.days / 30), "month")
    else:
        diff = delta.days / 365.0
        unit = '%0.1f years' % diff
    return "%s ago" % unit