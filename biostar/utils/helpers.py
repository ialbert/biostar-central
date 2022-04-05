from django.contrib.messages.storage import fallback
from django.test import RequestFactory, client
from django.conf import settings
import logging
import traceback
import html
import html2markdown
from datetime import datetime
from biostar import VERSION
import os
import uuid

logger = logging.getLogger('engine')


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def fake_request(url, data, user, method="POST", rmeta={}):
    "Make a fake request; defaults to POST."

    methods = {"POST": RequestFactory().post, "GET": RequestFactory().get,
               'PUT': RequestFactory().put}

    assert method in methods

    request = methods[method](url, data)

    # Mimic messaging system
    request.session = {}
    messages = fallback.FallbackStorage(request=request)
    request._messages = messages

    # Update the META info
    request.META.update(rmeta)

    request.user = user

    return request


def get_ip(request):
    """
    Attempts to extract the IP number from the HTTP request headers.
    """
    # lower the
    key = settings.IP_HEADER_KEY
    meta = request.META

    # Lowercase keys
    simple_meta = {k.lower(): v for k, v in request.META.items()}

    ip = meta.get(key, simple_meta.get(key, '0.0.0.0'))

    return ip


def htmltomarkdown(text):
    """
    Safely convert html to markdown
    """

    try:
        content = html2markdown.convert(text)
    except Exception as exc:
        logger.warning(f"error={exc};text={text[:100]}")
        # Return escaped text
        content = html.escape(text)

    return content


def ip_triplet(request):
    """
    Attempt to extract first three number from ip adress.
    """
    oip = get_ip(request=request)
    ips = oip.split(".")[:-1]
    ip = ".".join(ips)
    return ip
