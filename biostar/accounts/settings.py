from biostar.settings import *
from biostar.emailer.settings import *
# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

# Django debug flag.
DEBUG = True

# Should the site allow signup.
ALLOW_SIGNUP = False

ADMINS = [
    ("Admin User", "admin@localhost")
]

# The address the error emails are coming from.
SERVER_EMAIL = "admin@localhost"

# Whitelist of Ip addresses to not rate limit.
WHITELIST_IP = []

PAGEDOWN_APP = ['pagedown.apps.PagedownConfig']

PAGEDOWN_IMAGE_UPLOAD_ENABLED = True

# Allows admins to log in as any user.
DEBUG_USERS = True

# Maximum size per image uploaded, in mb.
MAX_IMAGE_SIZE_MB = 2

# Maximum number of images allowed.
MAX_IMAGES = 100

# User above this score do not get a reCAPTCHA
RECAPTCHA_THRESHOLD_USER_SCORE = 1

# The password for admin users. Must be changed in production.
DEFAULT_ADMIN_PASSWORD = "admin@localhost"

# Shortcut to first admin information.
ADMIN_NAME, ADMIN_EMAIL = ADMINS[0]

# The default sender name on emails.
DEFAULT_FROM_EMAIL = f"{ADMIN_NAME} <{ADMIN_EMAIL}>"

# User score threshold to be considered low reputation.
LOW_REP_THRESHOLD = 0

# Users below this threshold are considered to have recently joined.
RECENTLY_JOINED_DAYS = 30

# In MB
MAX_UPLOAD_SIZE = 10

# Trusted users upload limit in MB.
TRUSTED_UPLOAD_SIZE = 500

# Admin users upload limit in MB
ADMIN_UPLOAD_SIZE = 1000

MESSAGES_PER_PAGE = 5


# Additional middleware.
MIDDLEWARE += [
    #'biostar.accounts.middleware.limiter',
]


SIGNUP_RATE = '50/h'
PASSWORD_RESET_RATE = '50/h'

# Set RECAPTCH keys here.
RECAPTCHA_PUBLIC_KEY = ""
RECAPTCHA_PRIVATE_KEY = ""

# Django allauth settings.
SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True

# Key used to set ratelimitter.
# https://django-ratelimit.readthedocs.io/en/stable/security.html
# Must be set to the correct header
IP_HEADER_KEY = 'REMOTE_ADDR'
RATELIMIT_KEY = f"header:{IP_HEADER_KEY}"

# Rate to limit ( set to high value to disable ).
RATELIMIT_RATE = '5000/h'


# Other settings
ACCOUNT_AUTHENTICATION_METHOD = "email"
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_UNIQUE_EMAIL = True
ACCOUNT_USERNAME_REQUIRED = False
ACCOUNT_EMAIL_VERIFICATION = "mandatory"
ACCOUNT_EMAIL_SUBJECT_PREFIX = "[biostar] "
ACCOUNT_PASSWORD_MIN_LENGHT = 6
ACCOUNT_USER_MODEL_USERNAME_FIELD = None
ACCOUNT_USER_MODEL_EMAIL_FIELD = "email"

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

# Should the server look up locations in a task.
LOCATION_LOOKUP = True

INSTALLED_APPS = DEFAULT_APPS + ACCOUNTS_APPS + EMAILER_APP + PAGEDOWN_APP

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

    # ("Google", "A", "B"),
    # ("GitHub", "A", "B")
]

GOOGLE_TRACKER = ""
