from biostar.settings import *

import logging

logger = logging.getLogger("engine")

WSGI_APPLICATION = 'conf.test.test_wsgi.application'

ALLOWED_HOSTS += ['test.metabarcode.com']

HTTP_PORT = ''
PROTOCOL = 'https'

try:
    from .test_secrets import *
except ImportError as exc:
    logger.error(f"Module with local secrets could not be imported: {exc}")
