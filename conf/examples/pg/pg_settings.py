from biostar.forum.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ACCOUNTS_APPS + EMAILER_APP

DEBUG = True

WSGI_APPLICATION = 'conf.examples.pg.forum_wsgi.application'

DATABASE_NAME = os.environ.setdefault("DATABASE_NAME", "biostar-pgtest")

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

# Show debug toolbar
DEBUG_TOOLBAR = True

# Enable debug toolbar specific functions
if DEBUG_TOOLBAR:
    FORUM_APPS.append('debug_toolbar')
    MIDDLEWARE.append('debug_toolbar.middleware.DebugToolbarMiddleware')


try:
    from .postgres_secrets import *
except ImportError as exc:
    print("No postgres_secrets module could be imported")
