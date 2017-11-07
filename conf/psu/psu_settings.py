from biostar.settings import *

DEBUG = True

WSGI_APPLICATION = 'conf.psu.psu_wsgi.application'

SITE_ID = 1
SITE_DOMAIN = "psu.bioinformatics.recipes"
SITE_NAME = "PSU Bioinformatics Recipes"

ALLOWED_HOSTS = [SITE_DOMAIN]

SITE_HEADER = '<i class="barcode icon"></i> PSU Bioinformatics Recipes'

try:
    from .psu_secrets import *
except ImportError as exc:
    print("No recipes_secrets module could be imported")
