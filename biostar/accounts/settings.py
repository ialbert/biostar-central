from biostar.settings import *

import os

# Apply the logger settings.
from biostar.logconf import LOGGING

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

# Django debug flag.
DEBUG = True

# Override compression if needed.
# COMPRESS_ENABLED = True

# Change this in production!
SECRET_KEY = 'secret-key'

API_KEY = "api-key"

# Used internally by django to tested browser cookies
TEST_COOKIE_NAME = "test-cookie"

# Private key used to validate external logins
LOGIN_PRIVATE_KEY = "private-key"

# The password for admin users. Must be changed in production.
DEFAULT_ADMIN_PASSWORD = "admin@localhost"

RECAPTCHA_PUBLIC_KEY = ""
RECAPTCHA_PRIVATE_KEY = ""

SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True

# Set the home page to the engine or forum
INTERNAL_IPS = ['127.0.0.1']

# Admin users will be created automatically with DEFAULT_ADMIN_PASSWORD.
ADMINS = [
    ("Admin User", "admin@localhost")
]

# The default sender name on emails.
DEFAULT_FROM_EMAIL = f"{ADMINS[0][0]} <{ADMINS[0][1]}>"

# These must be set remote hosts.
SITE_ID = 1
SITE_DOMAIN = "localhost"
SITE_NAME = "Biostar Engine"

# Deployment specific parameters.
PROTOCOL = "http"
HTTP_PORT = ':8000'


BASE_URL = f"{PROTOCOL}://{SITE_DOMAIN}{HTTP_PORT}"


# Should the site allow signup.
ALLOW_SIGNUP = False

# Allow users to toggle their status moderator
ALLOW_SELF_MODERATE = False


LOGIN_REDIRECT_URL = "/project/list/private"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

# Helper function for building absolute paths.
def join(*args):
    return os.path.abspath(os.path.join(*args))


# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(join(__file__))

SOCIALACCOUNT_ADAPTER = "biostar.accounts.adapter.SocialAccountAdapter"

INSTALLED_APPS += [


    # Accounts configuration.
    'biostar.accounts.apps.AccountsConfig',

    # Allauth templates come last.
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
    'allauth.socialaccount.providers.github',
]


ROOT_URLCONF = 'biostar.accounts.urls'


AUTHENTICATION_BACKENDS += [
    "allauth.account.auth_backends.AuthenticationBackend",
]


# SENDFILE_BACKEND = "sendfile.backends.nginx"

SENDFILE_BACKEND = "sendfile.backends.development"

MEDIA_URL = '/media/'
STATIC_URL = '/static/'
STATIC_ROOT = join(BASE_DIR, '..', 'export', 'static')
STATICFILES_DIRS = [
    join(BASE_DIR, "engine", "static"),
    join(BASE_DIR, "forum", "static"),
]

STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    'compressor.finders.CompressorFinder',
)

# Apply default logger setting.
LOGGER_NAME = "engine"

# Use django-mailer to store emails in the database.
# EMAIL_BACKEND = "mailer.backend.DbBackend"

# The email delivery engine.
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'

# List of social login clients tuples.
# ( name, client_id, secret )

# Default clients redirect to localhost.
# Default clients may not be operational. See the
# django allauth documentataion on how to set them up.

#
# Callback example settings:
#
# http://localhost:8000/accounts/social/google/login/callback/
# http://localhost:8000/accounts/social/github/login/callback/
#
SOCIAL_CLIENTS = [

    ("Google", "547073349197-ri3eku9fdpi1ble7eoc4amlsrh7m2oiv.apps.googleusercontent.com", "DR1-zMqOLTqRGvhSDb5rQBMg"),
    ("GitHub", "d8493ce8967ea5abbd73", "04ff5043ebce2317d3c26bfe90fbc5e67fa38d05")

]
