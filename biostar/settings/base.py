# Django settings for biostar project.
import os, sys, re
from django.core.exceptions import ImproperlyConfigured
from .logger import LOGGING
from .social import *

# turn off debug mode on deployed servers
DEBUG = True

# template debug mode
TEMPLATE_DEBUG = DEBUG

def get_env(name):
    """Get the environment variable or return exception"""
    try:
        return os.environ[name]
    except KeyError:
        msg = "*** Required environment variable %s not set." % name
        raise ImproperlyConfigured(msg)

def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))

# displays debug comments when the server is run from this IP
INTERNAL_IPS = ('127.0.0.1', )

# the directory that this file is located in
__CURR_DIR = abspath(os.path.dirname(__file__))

# set location relative to the current file directory
HOME_DIR = abspath(__CURR_DIR, '..', '..')
DATABASE_DIR = abspath(HOME_DIR, 'data')
DATABASE_NAME = abspath(DATABASE_DIR, 'biostar2.db')
STATIC_DIR = abspath(HOME_DIR, 'biostar', 'static')
BIOSTAR_STATIC_ROOT = get_env("BIOSTAR_STATIC_ROOT")

# Must contains at least one (name, email) pair
ADMINS = (
    ('Istvan Albert', 'istvan.albert@gmail.com'),
)

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
ALLOWED_HOSTS = []

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
    # Put strings here, like "/home/html/static" or "C:/www/django/static".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
    STATIC_DIR,
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
#    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

# Make this unique, and don't share it with anybody.
# on deployed servers make this unique, and don't share it with anybody.
SECRET_KEY = get_env("SECRET_KEY")

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'biostar.urls'

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'biostar.wsgi.application'

TEMPLATE_DIRS = (
    # Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

LOGIN_REDIRECT_URL = "/"

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',

    # The javascript and CSS asset manager.
    'compressor',

    # Biostar specific apps.
    'biostar.apps.util',
    'biostar.apps.accounts',
    'biostar.apps.main',

    # Enabling the admin and its documentation
    'django.contrib.admin',
    'django.contrib.messages',
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    #'allauth.socialaccount.providers.github',
    'allauth.socialaccount.providers.google',
    #'allauth.socialaccount.providers.twitter',
    #'allauth.socialaccount.providers.facebook',
)

TEMPLATE_CONTEXT_PROCESSORS = (
    "django.core.context_processors.request",
    "django.contrib.auth.context_processors.auth",
    "allauth.account.context_processors.account",
    "allauth.socialaccount.context_processors.socialaccount",
)


# Session specific settings.
SESSION_SERIALIZER = 'django.contrib.sessions.serializers.JSONSerializer'

# Use a mock email backend for development.
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'