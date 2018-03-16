import logging
import os

logger = logging.getLogger("engine")

from django.core.wsgi import get_wsgi_application

# Override the DJANGO SETTINGS MODULE.
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "conf.natay.natay_settings")

application = get_wsgi_application()
