from biostar.forum.settings import *

import os
TRANSFER_APP = ['biostar.transfer']

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + PLANET_APPS + ACCOUNTS_APPS + EMAILER_APP + TRANSFER_APP

DEBUG = True

INIT_PLANET = False

WSGI_APPLICATION = 'conf.site.site_wsgi'

# The source database containing data from the old version of Biostar.
OLD_DATABASE = os.environ.setdefault("OLD_DATABASE", "old_biostar_db")

# The new database where the data will be copied into.
NEW_DATABASE = os.environ.setdefault("NEW_DATABASE", "database.db")
POSTGRES_HOST = os.environ.setdefault("POSTGRES_HOST", "")

print(f'NEW_DATABASE={NEW_DATABASE}, OLD_DATABASE={OLD_DATABASE}')

DATABASES = {

    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': NEW_DATABASE,
        'USER': '',
        'PASSWORD': '',
        'HOST': POSTGRES_HOST,
        'PORT': '',
    },

    'biostar2': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': OLD_DATABASE,
        'USER': '',
        'PASSWORD': '',
        'HOST': POSTGRES_HOST,
        'PORT': '',
        'TEST': {
            'MIRROR': 'default',
        }
    },
}

try:
    from conf.run.site_secrets import *
except ImportError as exc:
    print("No postgres_secrets module could be imported")
