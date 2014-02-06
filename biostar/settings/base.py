# -*- coding: utf8 -*-
#
# Django settings for biostar project.
#
import os
from django.core.exceptions import ImproperlyConfigured
from .logger import LOGGING

# Turn off debug mode on deployed servers.
DEBUG = True

# Template debug mode.
TEMPLATE_DEBUG = DEBUG


def get_env(name, func=None):
    """Get the environment variable or return exception"""
    try:
        if func:
            return func(os.environ[name])
        else:
            return unicode(os.environ[name], encoding="utf-8")
    except KeyError:
        msg = "*** Required environment variable %s not set." % name
        raise ImproperlyConfigured(msg)


def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))

# Displays debug comments when the server is run from this IP.
INTERNAL_IPS = ('127.0.0.1', )

# The directory that this file is located in.
__CURR_DIR = abspath(os.path.dirname(__file__))

# Set location relative to the current file directory.
HOME_DIR = get_env("BIOSTAR_HOME")
DATABASE_DIR = abspath(HOME_DIR, 'data')
DATABASE_NAME = abspath(DATABASE_DIR, 'biostar2.db')
STATIC_DIR = abspath(HOME_DIR, 'biostar', 'static')
BIOSTAR_STATIC_ROOT = get_env("BIOSTAR_STATIC_ROOT")
TEMPLATE_DIR = abspath(__CURR_DIR, '..', 'server', 'templates')

# Needs to point to the directory that contains the
# html files that are stored in the flatpages about, faq, help, policy etc.
FLATPAGE_IMPORT_DIR = abspath(TEMPLATE_DIR, 'flatpages')

# Default search index location.
DATA_DIR = abspath(__CURR_DIR, '..', '..', 'data')
WHOOSH_INDEX = abspath(DATA_DIR, "whoosh_index")

# These settings create an admin user.
# The default password is the SECRET_KEY.
ADMIN_NAME = get_env("BIOSTAR_ADMIN_NAME")
ADMIN_EMAIL = get_env("BIOSTAR_ADMIN_EMAIL")

ADMINS = (
    (ADMIN_NAME, ADMIN_EMAIL),
)

# Get the secret key from the environment.
SECRET_KEY = get_env("SECRET_KEY")

MANAGERS = ADMINS

DATABASES = {
    'default': {
        # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': DATABASE_NAME,
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    }
}

# admin site may fail if this setting is active
TEMPLATE_STRING_IF_INVALID = "*** MISSING ***"

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ["localhost"]



# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = 'America/Chicago'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

# These parameters will be inserted into the database automatically.
SITE_ID = 1
SITE_NAME = "localhost"
SITE_DOMAIN = "localhost"

DEFAULT_FROM_EMAIL = get_env("DEFAULT_FROM_EMAIL")

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale.
USE_L10N = True

# If you set this to False, Django will not use timezone-aware datetimes.
USE_TZ = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/var/www/example.com/media/"
MEDIA_ROOT = ''

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://example.com/media/", "http://media.example.com/"
MEDIA_URL = ''

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/var/www/example.com/static/"
STATIC_ROOT = BIOSTAR_STATIC_ROOT

# URL prefix for static files.
# Example: "http://example.com/static/", "http://static.example.com/"
STATIC_URL = '/static/'

# Additional locations of static files
STATICFILES_DIRS = (
    # Use absolute paths, not relative paths.
    STATIC_DIR,
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    'compressor.finders.CompressorFinder',
)

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'django.contrib.flatpages.middleware.FlatpageFallbackMiddleware',
    'biostar.server.middleware.Visit',
)

ROOT_URLCONF = 'biostar.urls'

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'biostar.wsgi.application'

TEMPLATE_DIRS = (
    TEMPLATE_DIR,
    # Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

LOGIN_REDIRECT_URL = "/"

MESSAGE_TAGS = {
    10: 'alert-info', 20: 'alert-info',
    25: 'alert-success', 30: 'alert-warning', 40: 'alert-error',
}

INSTALLED_APPS = [
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',

    # The javascript and CSS asset manager.
    'compressor',

    # Enabling the admin and its documentation.
    'django.contrib.sites',
    'django.contrib.admin',
    'django.contrib.messages',
    'django.contrib.humanize',
    'django.contrib.flatpages',

    # Social login handlers.
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
    'allauth.socialaccount.providers.twitter',
    'allauth.socialaccount.providers.facebook',
    'allauth.socialaccount.providers.persona',
    'allauth.socialaccount.providers.github',
    #'allauth.socialaccount.providers.linkedin',
    # 'allauth.socialaccount.providers.weibo',

    # External apps.
    'taggit',
    'haystack',
    'crispy_forms',

    # Biostar specific apps.
    'biostar.apps.util',
    'biostar.apps.posts',
    'biostar.apps.messages',
    'biostar.apps.badges',
    'biostar.apps.users',

    # The main Biostar server.
    'biostar.server',
]

CRISPY_TEMPLATE_PACK = 'bootstrap3'

AUTH_USER_MODEL = 'users.User'

DEBUG_TOOLBAR_PATCH_SETTINGS = False

# Default search is provided via Whoosh

HAYSTACK_CONNECTIONS = {
    'default': {
        'ENGINE': 'haystack.backends.whoosh_backend.WhooshEngine',
        'PATH': WHOOSH_INDEX,
    },
}

TEMPLATE_CONTEXT_PROCESSORS = (
    # Django specific context processors.
    "django.core.context_processors.debug",
    "django.core.context_processors.static",
    "django.core.context_processors.request",
    "django.contrib.auth.context_processors.auth",
    "django.contrib.messages.context_processors.messages",

    # Social authorization specific context.
    "allauth.account.context_processors.account",
    "allauth.socialaccount.context_processors.socialaccount",

    # Biostar specific context.
    'biostar.server.context.shortcuts',
)

AUTHENTICATION_BACKENDS = (
    'django.contrib.auth.backends.ModelBackend',
    "allauth.account.auth_backends.AuthenticationBackend",
)

ACCOUNT_CONFIRM_EMAIL_ON_GET = True

# Customize this to match the providers listed in the APPs
SOCIALACCOUNT_PROVIDERS = {

    'facebook': {
        'SCOPE': ['email'],
        'AUTH_PARAMS': {'auth_type': 'reauthenticate'},
        'METHOD': 'oauth2',
        'LOCALE_FUNC': lambda x: 'en_US',
        'PROVIDER_KEY': get_env("FACEBOOK_PROVIDER_KEY"),
        'PROVIDER_SECRET_KEY': get_env("FACEBOOK_PROVIDER_SECRET_KEY"),
    },

    'twitter': {
        'SCOPE': ['email'],
        'PROVIDER_KEY': get_env("TWITTER_PROVIDER_KEY"),
        'PROVIDER_SECRET_KEY': get_env("TWITTER_PROVIDER_SECRET_KEY"),
    },

    'persona': {
        'REQUEST_PARAMETERS': {'siteName': 'Biostar'}
    },

    'github': {
        'PROVIDER_KEY': get_env("GITHUB_PROVIDER_KEY"),
        'PROVIDER_SECRET_KEY': get_env("GITHUB_PROVIDER_SECRET_KEY"),
    },

    #'linkedin': {
    #    'SCOPE': ['r_emailaddress'],
    #    'PROFILE_FIELDS' : [ 'id', 'first-name', 'last-name', 'email-address'],
    #},

    #'weibo': {
    #    'PROVIDER_KEY': get_env("WEIBO_PROVIDER_KEY"),
    #    'PROVIDER_SECRET_KEY': get_env("WEIBO_PROVIDER_SECRET_KEY"),
    #},

    'google': {
        'SCOPE': ['email', 'https://www.googleapis.com/auth/userinfo.profile'],
        'AUTH_PARAMS': {'access_type': 'online'},
        'PROVIDER_KEY': get_env("GOOGLE_PROVIDER_KEY"),
        'PROVIDER_SECRET_KEY': get_env("GOOGLE_PROVIDER_SECRET_KEY"),
    },
}

# The google id will injected as a template variable.
GOOGLE_TRACKER = "foobar"

# The default CSS file to load.
SITE_STYLE_CSS = "biostar.style.less"

# Django precompressor settings.
COMPRESS_PRECOMPILERS = (
    ('text/coffeescript', 'coffee --compile --stdio'),
    ('text/less', 'lessc {infile} {outfile}'),
)

COMPRESS_OFFLINE_CONTEXT = {
    'STATIC_URL': STATIC_URL,
    'SITE_STYLE_CSS': SITE_STYLE_CSS,
}

# The cache mechanism is deployment dependent. Override it externally.
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.dummy.DummyCache' if DEBUG else 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'unique-snowflake'
    }
}

# This will form the category navigation bar.
# Will be treated as tags on posts.
# These should be the most frequent tags on the site.
DEFAULT_TOPICS = [
    "Forum", "Assembly", "RNA-Seq", "ChIP-Seq", "SNP-Calling", "Galaxy", "Jobs", "Planet",
]

# The number of posts to show per page.
PAGINATE_BY = 25

# How many posts to show per page.
POSTS_PER_PAGE = 15

ACCOUNT_AUTHENTICATION_METHOD = "email"
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_UNIQUE_EMAIL = True
ACCOUNT_USERNAME_REQUIRED = False
ACCOUNT_EMAIL_VERIFICATION = "optional"
ACCOUNT_EMAIL_SUBJECT_PREFIX = "[biostar] "
ACCOUNT_PASSWORD_MIN_LENGHT = 6
ACCOUNT_USER_MODEL_USERNAME_FIELD = None
ACCOUNT_USER_MODEL_EMAIL_FIELD = "email"
ACCOUNT_DEFAULT_HTTP_PROTOCOL = "http"
ACCOUNT_LOGOUT_ON_GET = True

# Session specific settings.
SESSION_SERIALIZER = 'django.contrib.sessions.serializers.JSONSerializer'

# Use a mock email backend for development.
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'

# On deployed servers the following must be set.
EMAIL_HOST = get_env("EMAIL_HOST")
EMAIL_PORT = get_env("EMAIL_PORT", func=int)
EMAIL_HOST_USER = get_env("EMAIL_HOST_USER")
EMAIL_HOST_PASSWORD = get_env("EMAIL_HOST_PASSWORD")