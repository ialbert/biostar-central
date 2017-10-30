from engine.settings import *
from .test_secrets import *

# This file will contain passwords that should not be in the repo.

WSGI_APPLICATION = 'conf.test.test_wsgi.application'

ALLOWED_HOSTS += ['test.metabarcode.com']
