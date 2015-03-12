from __future__ import absolute_import

VERSION = "3.0.0-alpha"

# Based on Celery recommendations.
from .celery import app as celery_app