from biostar.forum.settings import *
import os


DEBUG = True

CUSTOM_THEME = os.path.abspath(os.path.join(BASE_DIR, 'themes', 'bioconductor'))

HTTP_PROTOCOL = "http"

# Full file path to tags.
# Each line is a list of comma separated tags.
TAGS_OPTIONS_FILE = os.path.join(CUSTOM_THEME, 'tags', 'tags.txt')

REMAP_TYPE_DISPLAY = True
# Post types displayed when creating, empty list displays all types.
ALLOWED_POST_TYPES = ["Question", "Job", "Tutorial", "News"]


STATICFILES_DIRS = [os.path.join(CUSTOM_THEME, 'static')]

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
        'APP_DIRS':True,
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
WSGI_APPLICATION = "themes.bioconductor.conf.wsgi.application"

ADD_THREAD_USERS = False
GRAVATAR_ICON = 'identicon'

# These parameters will be inserted into the database automatically.
SITE_NAME = "Bioconductor Support Forum"

SITE_DOMAIN = "support.bioconductor.org"

# What domain will handle the replies.
EMAIL_REPLY_PATTERN = "reply+%s+code@bioconductor.org"

# The format of the email that is sent
EMAIL_FROM_PATTERN = u'''"%s [bioc]" <%s>'''

# The subject of the reply goes here
EMAIL_REPLY_SUBJECT = u"[bioc] %s"

SEARCH_LIMIT = 60

# On deployed servers the following must be set.
EMAIL_HOST = ""
EMAIL_PORT = ""
EMAIL_HOST_USER = ""
EMAIL_HOST_PASSWORD = ""
EMAIL_USE_TLS = False

AWS_ACCESS_KEY_ID = ''
AWS_SECRET_ACCESS_KEY = ''

EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
#EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'

#print(INDEX_DIR, INDEX_NAME, "INDEX")

# Tries to load up secret settings from a predetermined module
try:
    from themes.bioconductor.conf.secrets import *
    print(f"Loaded secrets from: themes.bioconductor.conf.run.secrets")
except Exception as exc:
    print(f"Secrets module not imported: {exc}")
