from .base import *

# add debugging middleware

MIDDLEWARE_CLASSES = list(MIDDLEWARE_CLASSES)
MIDDLEWARE_CLASSES.append('debug_toolbar.middleware.DebugToolbarMiddleware')

INSTALLED_APPS = list(INSTALLED_APPS)
INSTALLED_APPS.append('debug_toolbar')

