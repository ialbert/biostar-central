from biostar.settings import *

DEBUG = False

SITE_ID = 1
SITE_DOMAIN = "www.bioinformatics.recipes"
SITE_NAME = "Bioinformatics Recipes"

HTTP_PORT = ''
PROTOCOL = 'https'

ALLOWED_HOSTS = [SITE_DOMAIN]

WSGI_APPLICATION = 'conf.main.main_wsgi.application'

SITE_HEADER = '<i class="barcode icon"></i>Bioinformatics Recipes'

SECRET_KEY = '-replace-this-JksaLL8493^kasdaKKLSkjawalkj'

try:
    from .main_secrets import *
except ImportError as exc:
    print("No main_secrets module could be imported")
