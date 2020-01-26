from biostar.forum.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + PLANET_APPS + ACCOUNTS_APPS + EMAILER_APP

DEBUG = False

# Show debug toolbar
DEBUG_TOOLBAR = False

# Enable debug toolbar specific functions
if DEBUG_TOOLBAR:
    INSTALLED_APPS.extend([
        'debug_toolbar',
    ])
    MIDDLEWARE.append('debug_toolbar.middleware.DebugToolbarMiddleware')

# Add static serving capability
if DEBUG is False:
    print("Whitenoise static serve enabled (pip install whitenoise)")
    MIDDLEWARE.append('whitenoise.middleware.WhiteNoiseMiddleware')

POSTS_PER_PAGE = 50

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

# List additional directories from which templates should loaded from.
# First directory with the named template will load the template.
TEMPLATE_DIRS = [
    # join(BASE_DIR, "themes", "foo", "templates"),
]

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',

        'DIRS': TEMPLATE_DIRS,

        'APP_DIRS': True,

        'OPTIONS': {
            'string_if_invalid': "**MISSING**",
            'context_processors': [
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.template.context_processors.media',
                'django.contrib.messages.context_processors.messages',
                'biostar.context.main',
                'biostar.forum.context.forum',
            ],

            #'loaders': [
            #    ('django.template.loaders.cached.Loader', [
            #        'django.template.loaders.filesystem.Loader',
            #        'django.template.loaders.app_directories.Loader'
            #        ],)],

        },
    },
]

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'unique-snowflake',
    }
}

try:
    from .postgres_secrets import *
except ImportError as exc:
    print("No postgres_secrets module could be imported")

STATICFILES_FINDERS = [
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'compressor.finders.CompressorFinder'
]

STATICFILES_DIRS = [
    join(BASE_DIR, "biostar", "forum", "static"),
]


