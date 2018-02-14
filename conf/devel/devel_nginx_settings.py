from biostar.settings import *

DEBUG = False

WSGI_APPLICATION = 'conf.devel.devel_wsgi.application'

ALLOWED_HOSTS = ['localhost', 'www.lvh.me']

SITE_DOMAIN = "www.lvh.me"

HTTP_PORT = ":9000"

SENDFILE_BACKEND = 'sendfile.backends.nginx'

SENDFILE_ROOT = '/Users/ialbert/app/biostar-engine/export/media'
SENDFILE_URL = '/media'

try:
    from .devel_secrets import *
except ImportError as exc:
    print("No devel_secrets module could be imported")
