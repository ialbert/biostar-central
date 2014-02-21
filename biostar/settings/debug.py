# -*- coding: utf8 -*-
#
# Development settings
#
from .base import *

# add debugging middleware
MIDDLEWARE_CLASSES = list(MIDDLEWARE_CLASSES)
MIDDLEWARE_CLASSES.insert(0, 'debug_toolbar.middleware.DebugToolbarMiddleware')

# We load the debug toolbar as well
INSTALLED_APPS = list(INSTALLED_APPS)

# This needs to be added before the user models.
INSTALLED_APPS.append( "debug_toolbar")

def show_toolbar(request):
    return True

DEBUG_TOOLBAR_CONFIG ={
    'SHOW_TOOLBAR_CALLBACK': "biostar.apps.util.always_true",
}

# more email settings
EMAIL_USE_TLS = True

# Set up email backend.
EMAIL_BACKEND = 'biostar.apps.util.mailer.SSLEmailBackend'
EMAIL_HOST = get_env("EMAIL_HOST")
EMAIL_PORT = get_env("EMAIL_PORT", func=int)
EMAIL_HOST_USER = get_env("EMAIL_HOST_USER")
EMAIL_HOST_PASSWORD = get_env("EMAIL_HOST_PASSWORD")

XTEMPLATE_LOADERS = (
    (
        'django.template.loaders.cached.Loader', (
            'django.template.loaders.filesystem.Loader',
            'django.template.loaders.app_directories.Loader',
        )),
)