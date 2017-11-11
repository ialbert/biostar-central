from biostar.settings import *

WSGI_APPLICATION = 'conf.test.test_wsgi.application'

ALLOWED_HOSTS += ['test.metabarcode.com']

HTTP_PORT = ''
PROTOCOL = 'https'

try:
    from .test_secrets import *
except ImportError as exc:
    print("No test_secrets module could be imported")
