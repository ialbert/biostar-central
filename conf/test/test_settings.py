from biostar.engine.settings import *

try:
    from .test_secrets import *
except ImportError as err:
    print("No test_secrets found")

# This file will contain passwords that should not be in the repo.

WSGI_APPLICATION = 'conf.test.test_wsgi.application'

ALLOWED_HOSTS += ['test.metabarcode.com']

