# -*- coding: utf8 -*-
#
# Development settings
#
from .base import *

# These settings create an admin user. The default password is the SECRET_KEY.
ADMIN_NAME = u"Señor Admin"
ADMIN_EMAIL = u"foo@bar.com"

ADMINS = (
    (ADMIN_NAME, ADMIN_EMAIL),
)

# Get the secret key from the environment.
SECRET_KEY = get_env("SECRET_KEY")

# add debugging middleware
MIDDLEWARE_CLASSES = list(MIDDLEWARE_CLASSES)
MIDDLEWARE_CLASSES.insert(0, 'debug_toolbar.middleware.DebugToolbarMiddleware')

# We load the debug toolbar as well
INSTALLED_APPS = list(INSTALLED_APPS)

# This needs to be added before the user models.
INSTALLED_APPS.insert(0, "debug_toolbar")
