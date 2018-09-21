from biostar.settings import *

DEBUG = False

SITE_ID = 1
SITE_DOMAIN = "test.biostars.org"
SITE_NAME = "Test Biostar"

HTTP_PORT = ''
PROTOCOL = 'http'

ALLOWED_HOSTS = [SITE_DOMAIN]

WSGI_APPLICATION = 'conf.main.main_wsgi.application'

SITE_HEADER = '<i class="barcode icon"></i> Metagenomics Barcode Data Repository'

try:
    from .biostar_test_secrets import *
except ImportError as exc:
    print("No biostar_test_secrets module could be imported")
