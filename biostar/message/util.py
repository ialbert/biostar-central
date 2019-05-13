import re
import bleach
import logging
import uuid
from datetime import datetime


from django.template import loader
from django.utils.timezone import utc

logger = logging.getLogger(__name__)


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def render(name, **kwds):
    "Helper function to render a template"
    tmpl = loader.get_template(name)
    page = tmpl.render(kwds)
    return page


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]
