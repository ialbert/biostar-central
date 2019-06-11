from biostar.engine.settings import *
from biostar.forum.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ENGINE_APPS + ACCOUNTS_APPS + EMAILER_APP

ROOT_URLCONF = 'biostar.test.test_urls'

DEBUG = True


WSGI_APPLICATION = 'conf.examples.postgres.postgres_wsgi.application'

DATABASES = {

    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'engine.db',
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
