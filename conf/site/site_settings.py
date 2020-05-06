import os
import uuid

from biostar.settings import *

# from biostar.recipes.settings import *

from biostar.forum.settings import *

import logging
import platform

logger = logging.getLogger("biostar")

DEBUG = True

# Set your own secret key here.
SECRET_KEY = str(uuid.uuid4())

# Admin users will be created automatically with DEFAULT_ADMIN_PASSWORD.
ADMINS = [
    ("Admin User", "admin@localhost")
]

# Set the default adming password.
DEFAULT_ADMIN_PASSWORD = SECRET_KEY

# Attempts to detect hostname automatically.
# Override with your own domain.
SITE_DOMAIN = platform.node()

SITE_ID = 1

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
