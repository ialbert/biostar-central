# -*- coding: utf8 -*-
#
# Deployment settings
#
from biostar.settings.base import *

DEBUG = False
TEMPLATE_DEBUG = DEBUG

TEMPLATE_LOADERS = (
    (
        'django.template.loaders.cached.Loader', (
        'django.template.loaders.filesystem.Loader',
        'django.template.loaders.app_directories.Loader',
    )),
)

# more email settings
EMAIL_USE_TLS = True

# Set up email backend.
EMAIL_BACKEND = 'biostar.apps.util.mailer.SSLEmailBackend'
EMAIL_HOST = get_env("EMAIL_HOST")
EMAIL_PORT = get_env("EMAIL_PORT", func=int)
EMAIL_HOST_USER = get_env("EMAIL_HOST_USER")
EMAIL_HOST_PASSWORD = get_env("EMAIL_HOST_PASSWORD")

