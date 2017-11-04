from biostar.engine.settings import *

DEBUG = True

WSGI_APPLICATION = 'conf.recipes.recipes_wsgi.application'

SITE_ID = 1
SITE_DOMAIN = "www.bioinformatics.recipes"
SITE_NAME = "Bioinformatics Recipes"

ALLOWED_HOSTS = [ SITE_DOMAIN  ]

SITE_HEADER =  '<i class="barcode icon"></i> Bioinformatics Recipes'

try:
    from .recipes_secrets import *
except ImportError as exc:
    print("No recipes.recipes_secrets module could be found")
