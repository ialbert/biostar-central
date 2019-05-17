from biostar.engine.settings import *

from biostar.settings import *
from biostar.accounts.settings import *
from biostar.message.settings import *
from biostar.emailer.settings import *
from biostar.engine.settings import *
from biostar.forum.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ENGINE_APPS + MESSAGE_APPS + ACCOUNTS_APPS  + EMAILER_APP

ROOT_URLCONF = 'biostar.test.test_urls'

DEBUG = True

WSGI_APPLICATION = 'conf.postgres.postgres_wsgi.application'

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
