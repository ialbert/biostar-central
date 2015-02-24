"""
Handles the output formatting and conversions
"""
import bleach
from django.conf import settings
from markdown2 import markdown

# Tags that are allowed by default
ALLOWED_TAGS = bleach.ALLOWED_TAGS + settings.ALLOWED_TAGS
ALLOWED_STYLES = bleach.ALLOWED_STYLES + settings.ALLOWED_STYLES
ALLOWED_ATTRIBUTES = dict(bleach.ALLOWED_ATTRIBUTES)
ALLOWED_ATTRIBUTES.update(settings.ALLOWED_ATTRIBUTES)

# Moderators may use more tags and styles
TRUSTED_TAGS = ALLOWED_TAGS + settings.TRUSTED_TAGS
TRUSTED_STYLES = ALLOWED_STYLES + settings.TRUSTED_STYLES
TRUSTED_ATTRIBUTES = dict(ALLOWED_ATTRIBUTES)
#settings.TRUSTED_ATTRIBUTES

def clean(text):
    "Sanitize text with no other substitutions"
    result = bleach.clean(text, tags=ALLOWED_TAGS, attributes=ALLOWED_ATTRIBUTES, styles=ALLOWED_STYLES)
    return result

def generate(text):
    "Sanitize text and expand links to match content"

    # To avoid circular imports
    from biostar3.forum.models import User, Post

    callbacks = []

    html = bleach.clean(text, tags=ALLOWED_TAGS, attributes=ALLOWED_ATTRIBUTES, styles=ALLOWED_STYLES)

    html = markdown(html)

    html = bleach.linkify(html, callbacks=callbacks, skip_pre=True)

    return html