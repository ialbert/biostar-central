from biostar.forum.settings import *
import os
from biostar.forum.settings import BASE_DIR

DEBUG = True

CUSTOM_THEME = os.path.abspath(os.path.join(BASE_DIR, 'themes', 'bioconductor'))

ACCOUNT_EMAIL_SUBJECT_PREFIX = "[bioc] "
HTTP_PROTOCOL = "http"

LANGUAGE_DETECTION = ["en"]

# Full file path to tags.
# Each line is a list of comma separated tags.
TAGS_OPTIONS_FILE = os.path.join(CUSTOM_THEME, 'tags', 'tags.txt')

# Ensure at least one tag in file is included.
REQUIRED_TAGS = os.path.join(CUSTOM_THEME, 'tags', 'packageList.txt')

REQUIRED_TAGS_URL = 'http://bioconductor.org/packages/devel/BiocViews.html#___Software'

# Rate to limit
RATELIMIT_RATE = '200/d'

# Post types displayed when creating, empty list displays all types.
ALLOWED_POST_TYPES = ["Question", "Job", "Tutorial", "News"]

DOCS_ROOT = os.path.join(CUSTOM_THEME, 'docs')

STATICFILES_DIRS = [os.path.join(CUSTOM_THEME, 'static'), DOCS_ROOT] + STATICFILES_DIRS

if DEBUG:
    TEMPLATE_LOADERS = (
        'django.template.loaders.filesystem.Loader',
        'django.template.loaders.app_directories.Loader',
    )
else:
    TEMPLATE_LOADERS = (
        'django.template.loaders.cached.Loader',
        'django.template.loaders.filesystem.Loader',
        'django.template.loaders.app_directories.Loader',
    )

# Template specific settings.
TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(CUSTOM_THEME, 'templates')],
        'APP_DIRS': True,
        'OPTIONS': {
            'string_if_invalid': "**MISSING**",
            'context_processors': [
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.template.context_processors.media',
                'django.contrib.messages.context_processors.messages',
                'biostar.forum.context.forum',
            ],
        },
    },
]

DOCS_ROOT = os.path.join(CUSTOM_THEME, 'docs')

# How many recent objects to show in the feed.
VOTE_FEED_COUNT = 5
LOCATION_FEED_COUNT = 5
AWARDS_FEED_COUNT = 5
REPLIES_FEED_COUNT = 5

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'biostar.wsgi.application'

CLASSIFY_SPAM = False

ADD_THREAD_USERS = False
GRAVATAR_ICON = 'identicon'

# These parameters will be inserted into the database automatically.
SITE_NAME = "Bioconductor Support Forum"

SITE_DOMAIN = "127.0.0.1"

# What domain will handle the replies.
EMAIL_REPLY_PATTERN = "reply+%s+code@bioconductor.org"

STRICT_TAGS = False

# The format of the email that is sent
EMAIL_FROM_PATTERN = u'''"%s [bioc]" <%s>'''

# The subject of the reply goes here
EMAIL_REPLY_SUBJECT = u"[bioc] %s"

# The default no reply email.
DEFAULT_NOREPLY_EMAIL = "noreply@bioconductor.org"

SEARCH_LIMIT = 60

# On deployed servers the following must be set.
EMAIL_HOST = ""
EMAIL_PORT = ""
EMAIL_HOST_USER = ""
EMAIL_HOST_PASSWORD = ""
EMAIL_USE_TLS = False

AWS_ACCESS_KEY_ID = ''
AWS_SECRET_ACCESS_KEY = ''

DROPDOWN_TAGS = True
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'

# Tries to load up secret settings from a predetermined module
try:
    from conf.run.site_secrets import *

    print(f"Loaded secrets from: themes.bioconductor.conf.run.secrets")
except Exception as exc:
    print(f"Secrets module not imported: {exc}")
