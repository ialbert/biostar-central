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
