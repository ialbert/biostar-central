from biostar.recipes.settings import *

INSTALLED_APPS = DEFAULT_APPS + ENGINE_APPS + ACCOUNTS_APPS + EMAILER_APP

DEBUG = False

# Show debug toolbar
DEBUG_TOOLBAR = False

# Enable debug toolbar specific functions
if DEBUG_TOOLBAR:
    INSTALLED_APPS.extend([
        'debug_toolbar',
    ])
    MIDDLEWARE.append('debug_toolbar.middleware.DebugToolbarMiddleware')
    print("FOOOO BARRR")

# Add static serving capability
if DEBUG is False:
    print("Whitenoise static serve enabled (pip install whitenoise)")
    MIDDLEWARE.append('whitenoise.middleware.WhiteNoiseMiddleware')

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
