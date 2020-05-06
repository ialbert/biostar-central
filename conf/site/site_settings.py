import os

from biostar.settings import *

# from biostar.recipes.settings import *

from biostar.forum.settings import *

import logging
import platform

logger = logging.getLogger("biostar")

DEBUG = True

SITE_ID = 1

# Attempts to detect hostname automatically.
SITE_DOMAIN = platform.node()
SITE_NAME = "Biostar Central"
HTTP_PORT = ''
PROTOCOL = 'http'

ALLOWED_HOSTS = [SITE_DOMAIN]

DATABASE_NAME = "biostar-database"

DATABASES = {

    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': DATABASE_NAME,
        'USER': '',
        'PASSWORD': '',
        'HOST': '/var/run/postgresql/',
        'PORT': '',
    },
}

WSGI_APPLICATION = 'conf.run.site_wsgi.application'

try:
    # Attempts to load site secrets.
    from .site_secrets import *

    logger.info("Imported settings from .site_secrets")
except ImportError as exc:
    logger.warn(f"No secrets module could be imported: {exc}")
