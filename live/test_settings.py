from engine.settings import *

# This file will contain passwords that should not be in the repo.
from live.test_secrets import *

WSGI_APPLICATION = 'live.test_wsgi.application'

ALLOWED_HOSTS += [ 'test.metabarcode.com' ]

