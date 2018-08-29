from biostar.settings import *

#DEBUG = False

#SITE_ID = 1
#SITE_DOMAIN = "www.bioinformatics.recipes"
#SITE_NAME = "Bioinformatics Recipes"

#HTTP_PORT = ''
#PROTOCOL = 'https'

#ALLOWED_HOSTS = [SITE_DOMAIN]

WSGI_APPLICATION = 'conf.main.main_wsgi.application'

SITE_HEADER = '<i class="barcode icon"></i>Bioinformatics Recipes'

SECRET_KEY = "foo"

try:
    from .site_secrets import *
    print("Imported settings from '.site_secrets")
except ImportError as exc:
    print("No site_secrets module could be imported")
