from biostar.settings import *
from biostar.emailer.settings import *
# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

# Django debug flag.
DEBUG = True

# Should the site allow signup.
ALLOW_SIGNUP = False

# Private key used to validate external logins
LOGIN_PRIVATE_KEY = "private-key"

# In MB
MAX_UPLOAD_SIZE = 10

MESSAGES_PER_PAGE = 5
# The password for admin users. Must be changed in production.
DEFAULT_ADMIN_PASSWORD = "admin@localhost"

RECAPTCHA_PUBLIC_KEY = ""
RECAPTCHA_PRIVATE_KEY = ""

SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True


LOGIN_REDIRECT_URL = "/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

SOCIALACCOUNT_ADAPTER = "biostar.accounts.adapter.SocialAccountAdapter"

ACCOUNTS_APPS = [

    # Accounts configuration.
    'biostar.accounts.apps.AccountsConfig',

    # Allauth templates come last.
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
    'allauth.socialaccount.providers.github',
]

INSTALLED_APPS = DEFAULT_APPS + ACCOUNTS_APPS + EMAILER_APP

AUTHENTICATION_BACKENDS += ["allauth.account.auth_backends.AuthenticationBackend"]

ROOT_URLCONF = 'biostar.accounts.urls'

# List of social login clients tuples.
# ( name, client_id, secret )

# Default clients redirect to localhost.
# Default clients may not be operational. See the
# django allauth documentation on how to set them up.

#
# Callback example settings:
#
# http://localhost:8000/accounts/social/google/login/callback/
# http://localhost:8000/accounts/social/github/login/callback/
#
SOCIAL_CLIENTS = [
    ("Google", "A", "B"),
    ("GitHub", "A", "B")
]

try:
    from conf.secrets.defaults import *
    print (f"*** loaded 'conf.secrets.defaults'")
except ImportError as err:
    print (f"*** unable to import secrets: {err}")


