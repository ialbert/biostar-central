from biostar.settings import *

DEBUG = False

WSGI_APPLICATION = 'conf.natay.natay_wsgi.application'

SITE_ID = 1
SITE_DOMAIN = "natay.software"
SITE_NAME = "Software Demo"

HTTP_PORT = ''
PROTOCOL = 'https'

ALLOWED_HOSTS = [SITE_DOMAIN, "165.227.92.172"]

SITE_HEADER = '<img class="ui middle aligned small image" src="/static/images/shield.png"> <span>Demo the Biostar Engine</span>'


SENDFILE_BACKEND = 'sendfile.backends.nginx'

SENDFILE_ROOT = '/biostar-engine/export/media'
SENDFILE_URL = '/media'


try:
    from .psu_secrets import *
except ImportError as exc:
    print("No recipes_secrets module could be imported")
