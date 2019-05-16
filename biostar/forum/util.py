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

