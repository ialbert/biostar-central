from biostar.forum.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ACCOUNTS_APPS + EMAILER_APP

DEBUG = True

# Show debug toolbar
DEBUG_TOOLBAR = True

# Enable debug toolbar specific functions
if DEBUG_TOOLBAR:
    INSTALLED_APPS.append('debug_toolbar')
    MIDDLEWARE.append('debug_toolbar.middleware.DebugToolbarMiddleware')

WSGI_APPLICATION = 'conf.examples.pg.pg_wsgi.application'

DATABASE_NAME = os.getenv("DATABASE_NAME", "database.db")

DATABASES = {

    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': DATABASE_NAME,
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    },
}

try:
    from .postgres_secrets import *
except ImportError as exc:
    print("No postgres_secrets module could be imported")
