from biostar.settings import *

# Enable the right settings.

#from biostar.recipes.settings import *

from biostar.forum.settings import *

import logging

logger = logging.getLogger("biostar")

DEBUG = False

# Edit the site domain for specific cases.
SITE_DOMAIN = "www.lvh.me"

SITE_ID = 1

SITE_NAME = "Biostar Engine"
HTTP_PORT = ''
PROTOCOL = 'http'

ALLOWED_HOSTS = [SITE_DOMAIN]

WSGI_APPLICATION = 'conf.run.site_wsgi.application'

try:
    # Attempts to load site secrets.
    from .site_secrets import *

    logger.info("Imported settings from .site_secrets")
except ImportError as exc:
    logger.warn(f"No secrets module could be imported: {exc}")
