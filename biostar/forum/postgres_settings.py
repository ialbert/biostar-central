from biostar.settings import *

DEBUG = True

ONLY_FORUM_URLS = True
ENGINE_AS_ROOT = False

MENU_BAR = "widgets/forum_menubar.html"

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
