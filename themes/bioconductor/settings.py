from biostar.forum.settings import *
import socket
import os

#
# Override settings as needed
#


DEBUG = False

# Define Categories in top banner
START_CATEGORIES = ["Latest", "News", "Jobs", "Tutorials"]

NAVBAR_TAGS = []

END_CATEGORIES = []

# These are the tags that always show up in the tag recommendation dropdown.
POST_TAG_LIST = NAVBAR_TAGS + ["software error"]

# This will form the navbar
CATEGORIES = START_CATEGORIES + NAVBAR_TAGS + END_CATEGORIES

# This will appear as a top banner.
# It should point to a template that will be included.
TOP_BANNER = ""
# TOP_BANNER = "bioc_banner.html"

# Custom directory with bioconductor theme/
CUSTOM_THEME = os.path.abspath(os.path.join(BASE_DIR, '..', 'themes', 'bioconductor'))

STATICFILES_DIRS += [os.path.join(CUSTOM_THEME, 'static')]

# Template specific settings.
TEMPLATES[0]['TEMPLATE_DIRS'] += [os.path.join(CUSTOM_THEME, 'templates')]

# The site logo image on top navigation bar
SITE_LOGO = "bioconductor_logo_rgb_small.jpg"

# How many recent objects to show in the feed.
VOTE_FEED_COUNT = 5
LOCATION_FEED_COUNT = 5
AWARDS_FEED_COUNT = 5
REPLIES_FEED_COUNT = 5

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'biostar.wsgi.application'

# These parameters will be inserted into the database automatically.
SITE_NAME = "support.bioconductor.org"

# What domain will handle the replies.
EMAIL_REPLY_PATTERN = "reply+%s+code@bioconductor.org"

# The format of the email that is sent
EMAIL_FROM_PATTERN = u'''"%s [bioc]" <%s>'''

# The subject of the reply goes here
EMAIL_REPLY_SUBJECT = u"[bioc] %s"

# List of callables that know how to import templates from various sources.
if not DEBUG:
    TEMPLATES[0]['OPTIONS']['loaders'] = [
        (
            'django.template.loaders.cached.Loader',
            (
                'django.template.loaders.filesystem.Loader',
                'django.template.loaders.app_directories.Loader',
            )
        ),
    ]

# Other settings
ACCOUNT_AUTHENTICATION_METHOD = "email"
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_UNIQUE_EMAIL = True
ACCOUNT_USERNAME_REQUIRED = False
ACCOUNT_EMAIL_VERIFICATION = "mandatory"
ACCOUNT_EMAIL_SUBJECT_PREFIX = "[bioc] "
ACCOUNT_PASSWORD_MIN_LENGHT = 6
ACCOUNT_USER_MODEL_USERNAME_FIELD = None
ACCOUNT_USER_MODEL_EMAIL_FIELD = "email"
if ("ON_PRODUCTION" in os.environ and os.environ['ON_PRODUCTION'] == 'True'):
    ACCOUNT_DEFAULT_HTTP_PROTOCOL = "https"
else:
    ACCOUNT_DEFAULT_HTTP_PROTOCOL = "http"

# Google ReCaptcha No-Captcha settings
# When set the captcha forms will be active.
if "RECAPTCHA_PUBLIC_KEY" in os.environ:
    RECAPTCHA_PUBLIC_KEY = os.environ['RECAPTCHA_PUBLIC_KEY']
else:
    RECAPTCHA_PUBLIC_KEY = ""
if "RECAPTCHA_PRIVATE_KEY" in os.environ:
    RECAPTCHA_PRIVATE_KEY = os.environ['RECAPTCHA_PRIVATE_KEY']
else:
    RECAPTCHA_PRIVATE_KEY = ""


RECAPTCHA_USE_SSL = True  # Defaults to False
NOCAPTCHA = True


# Email
if socket.gethostname() in ["gamay", "habu"] or ("EMAIL_HOST" in os.environ.keys() and (get_env("EMAIL_HOST") == "mailcatcher") \
         or ("amazonaws.com" in get_env("EMAIL_HOST").lower())):
    EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'
else:
    EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'

# On deployed servers the following must be set.
EMAIL_HOST = get_env("EMAIL_HOST")
EMAIL_PORT = get_env("EMAIL_PORT", func=int)
EMAIL_HOST_USER = get_env("EMAIL_HOST_USER")
EMAIL_HOST_PASSWORD = get_env("EMAIL_HOST_PASSWORD")
EMAIL_USE_TLS = False

if ("EMAIL_USE_TLS" in os.environ.keys() and get_env("EMAIL_USE_TLS").lower() == "true"):

    EMAIL_USE_TLS=True

# Tries to load up secret settings from a predetermined module
try:
    from conf.run.secrets import *

    print(f"Loaded secrets from: conf.run.secrets")
except Exception as exc:
    print(f"Secrets module not imported: {exc}")
