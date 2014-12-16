from __future__ import absolute_import
import os

from django.core.exceptions import ImproperlyConfigured
from .logger import LOGGING


def abspath(*args):
    "Generates absolute paths."
    return os.path.abspath(os.path.join(*args))


def get_env(name, func=None):
    "Gets the values from environment variables. Applies the function if specified."
    try:
        if func:
            return func(os.environ[name])
        else:
            return unicode(os.environ[name], encoding="utf-8")
    except KeyError:
        msg = "*** Required environment variable %s not set. See conf/defaults.env." % name
        raise ImproperlyConfigured(msg)

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os

BASE_DIR = get_env('BIOSTAR_HOME')


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.7/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = get_env('SECRET_KEY')

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

TEMPLATE_DEBUG = True

ALLOWED_HOSTS = []


# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.sites',
    'django.contrib.sitemaps',
    'biostar3.forum',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'biostar3.urls'

WSGI_APPLICATION = 'biostar3.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.7/ref/settings/#databases
DATABASE_NAME = get_env('DATABASE_NAME')
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': DATABASE_NAME,
    }
}

AUTH_USER_MODEL = 'forum.User'

TEMPLATE_PATH = abspath(get_env('TEMPLATE_PATH'))
DEFAULT_PATH = abspath(BASE_DIR, "biostar3", "themes", "default")
TEMPLATE_DIRS = (
   TEMPLATE_PATH,
   DEFAULT_PATH
)

# Internationalization
# https://docs.djangoproject.com/en/1.7/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.7/howto/static-files/

STATIC_URL = '/static/'
STATIC_ROOT = get_env('STATIC_ROOT')
STATICFILES_DIRS = (
    abspath(TEMPLATE_PATH, "static"),
)