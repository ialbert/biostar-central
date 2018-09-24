import re
import bleach
import logging
import requests
import mistune
import uuid
from datetime import datetime
#from html5lib.tokenizer import HTMLTokenizer


from django.template import loader, Context
from django.utils.timezone import utc

logger = logging.getLogger(__name__)


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
    text = bleach.clean(text, tags=[], attributes=[], styles={}, strip=True)
    return text
