from biostar.engine.settings import *
from biostar.forum.settings import *
import os
TRANSFER_APP = ['biostar.transfer']

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ENGINE_APPS + ACCOUNTS_APPS + EMAILER_APP + TRANSFER_APP

ROOT_URLCONF = 'biostar.test.test_urls'

DEBUG = True

EMAIL_BACKEND = "django.core.mail.backends.console.EmailBackend"

WSGI_APPLICATION = 'conf.examples.postgres.postgres_wsgi.application'

DATABASE_NAME = os.environ.setdefault("DATABASE_NAME", "database.db")

DATABASES = {

    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': DATABASE_NAME,
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    },

    'biostar2': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'biostar.db',
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
        'TEST': {
            'MIRROR': 'default',
        }
    },
}

try:
    from .postgres_secrets import *
except ImportError as exc:
    print("No postgres_secrets module could be imported")
