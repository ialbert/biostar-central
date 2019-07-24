from biostar.settings import *

import logging

logger = logging.getLogger("biostar")

DEBUG = False

SITE_ID = 1
SITE_DOMAIN = "www.bioinformatics.recipes"
SITE_NAME = "Bioinformatics Recipes"

HTTP_PORT = ''
PROTOCOL = 'https'

ALLOWED_HOSTS = [SITE_DOMAIN]

WSGI_APPLICATION = 'conf.run.site_wsgi.application'

try:
    # Attempts to load site secrets.
    from .site_secrets import *
    logger.info("Imported settings from .site_secrets")
except ImportError as exc:
    logger.warn(f"No secrets module could be imported: {exc}")
