"""
Constants that may be used in multiple packages
"""
import bleach

ALLOWED_TAGS = bleach.ALLOWED_TAGS + "p div br code pre".split()
ALLOWED_STYLES = bleach.ALLOWED_STYLES
ALLOWED_ATTRIBUTES = bleach.ALLOWED_ATTRIBUTES

LOCAL_MESSAGE, EMAIL_MESSAGE = range(2)
MESSAGING_TYPE_CHOICES = ((LOCAL_MESSAGE, "Local"), (EMAIL_MESSAGE, "Email"))