import os, logging
from django.core.wsgi import get_wsgi_application
from django.core.management import call_command

logger= logging.getLogger("biostar")

# Override the DJANGO SETTINGS MODULE

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "themes.bioconductor.settings")

dj = os.environ["DJANGO_SETTINGS_MODULE"]

print(f"*** DJANGO_SETTINGS_MODULE={dj}")
print('**** Mounting app ')
application = get_wsgi_application()
