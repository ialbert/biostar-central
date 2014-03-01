# -*- coding: utf8 -*-
#
# Deployment settings
#
from biostar.settings.base import *

DEBUG = False
TEMPLATE_DEBUG = DEBUG

# How far to look for posts for anonymous users.
COUNT_INTERVAL_WEEKS = 10

# How frequently do we update the counts for authenticated users.
SESSION_UPDATE_SECONDS = 10 * 60

# The top navigation has three parts: start, special tags, end
NAVBAR_START_CATEGORIES = [
    "Latest", "Unanswered",
]

# These should be the most frequent (or special) tags on the site.
NAVBAR_SPECIAL_TAGS = [
    "Assembly", "RNA-Seq", "ChIP-Seq", "SNP-Calling", "Galaxy",
]

NAVBAR_END_CATEGORIES = [
    "Job", "Planet", "Forum",
]

# This will form the navbar
CATEGORIES = NAVBAR_START_CATEGORIES + NAVBAR_SPECIAL_TAGS + NAVBAR_END_CATEGORIES

# Enable this if you have the lessc installed.
USE_COMPRESSOR = True

# The celery configuration file
CELERY_CONFIG = 'live.celeryconfig'

#BROKER_URL = 'redis://localhost:6379/0'
#CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'

ADMINS = (
    (ADMIN_NAME, ADMIN_EMAIL),
)

# These parameters will be inserted into the database automatically.
SITE_ID = 1
SITE_NAME = "My Site"

# The google id will injected as a template variable.
GOOGLE_TRACKER = "foobar"

# The default CSS file to load.
SITE_STYLE_CSS = "biostar.style.less"

# The site logo.
SITE_LOGO = "biostar.logo.png"

# Load the database name from the environment.
DATABASE_NAME = get_env("DATABASE_NAME")

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': DATABASE_NAME,
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    }
}


# Setting a cookie with email:signed_hash(email)
# will automatically create accounts
EXTERNAL_AUTH = [
    ("foo.bar.com", "ABC"),
]

# This will speed up the templating system.
TEMPLATE_LOADERS = (
    (
        'django.template.loaders.cached.Loader', (
            'django.template.loaders.filesystem.Loader',
            'django.template.loaders.app_directories.Loader',
        )),
)

# Background email sending.
EMAIL_USE_TLS = True

CELERY_EMAIL_BACKEND = 'biostar.mailer.SSLEmailBackend'

EMAIL_BACKEND = 'biostar.mailer.CeleryEmailBackend'
