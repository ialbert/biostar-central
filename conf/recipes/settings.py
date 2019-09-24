from biostar.settings import *

from biostar.recipes.settings import *
import logging

logger = logging.getLogger("engine")

SITE_DOMAIN = 'test.bioinformatics.recipes'

WSGI_APPLICATION = 'conf.recipes.wsgi.application'

ALLOWED_HOSTS += [SITE_DOMAIN]

HTTP_PORT = ''
PROTOCOL = 'https'

try:
    from .secrets import *
except ImportError as exc:
    logger.error(f"Module with local secrets could not be imported: {exc}")
