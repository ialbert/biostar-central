import os, logging
from django.core.wsgi import get_wsgi_application
from django.core.management import call_command

logger = logging.getLogger("biostar")

# Detect the active Django settings module.

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "conf.run.site_settings")

value = os.environ["DJANGO_SETTINGS_MODULE"]

print(f"*** DJANGO_SETTINGS_MODULE={value}")

application = get_wsgi_application()
