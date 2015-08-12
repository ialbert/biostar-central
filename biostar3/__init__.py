import sys

VERSION = "3.0.0-alpha"

# Based on Celery recommendations.
from .celery import app as celery_app


if sys.version_info < (3, 4):
    raise SystemExit("\n*** This version of Biostar requires Python 3.4 or greater.\n")