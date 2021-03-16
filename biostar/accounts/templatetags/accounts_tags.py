
import logging

from datetime import timedelta
from django import template
from django.conf import settings

from biostar.accounts.util import now
logger = logging.getLogger("engine")
register = template.Library()


@register.inclusion_tag('widgets/show_messages.html')
def show_messages(messages):
    """
    Renders the messages
    """
    return dict(messages=messages)


@register.simple_tag
def http_port():
    return f":{settings.HTTP_PORT}" if settings.HTTP_PORT else ""


@register.filter
def show_email(target, user=None):

    # Show the email to the same user.
    if target == user:
        return target.email

    try:
        head, tail = target.email.split("@")
        email = head[0] + "*" * len(tail) + "@" + tail
    except Exception as exc:
        return "*" * 10

    return email



