from __future__ import absolute_import
import os, string
from django.core.exceptions import ImproperlyConfigured

# Logging configuration.
from .logger import LOGGING

# This pulls in various site specific settings.
from .values import *

def abspath(*args):
    "Generates absolute paths."
    return os.path.abspath(os.path.join(*args))

def get_env(name, default=''):
    "Gets values from environment variables."
    value = os.environ.get(name, default)
    if not value:
        msg = "*** Required environment variable %s not set. See README.md " % name
        raise ImproperlyConfigured(msg)
    return value

TEST_RUNNER = 'django.test.runner.DiscoverRunner'

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
    'django.contrib.humanize',
    'compressor',
    'haystack',
    'taggit',
    'biostar3.forum',

    # Authentication apps
    'allauth',
    'allauth.account',
    'allauth.socialaccount',

    # Available providers:
    'allauth.socialaccount.providers.persona',
    'allauth.socialaccount.providers.github',
    'allauth.socialaccount.providers.google',

    # 'allauth.socialaccount.providers.amazon',
    # 'allauth.socialaccount.providers.angellist',
    # 'allauth.socialaccount.providers.bitbucket',
    # 'allauth.socialaccount.providers.bitly',
    # 'allauth.socialaccount.providers.coinbase',
    # 'allauth.socialaccount.providers.dropbox',
    # 'allauth.socialaccount.providers.dropbox_oauth2',
    # 'allauth.socialaccount.providers.facebook',
    # 'allauth.socialaccount.providers.flickr',
    # 'allauth.socialaccount.providers.feedly',
    # 'allauth.socialaccount.providers.fxa',
    # 'allauth.socialaccount.providers.hubic',
    # 'allauth.socialaccount.providers.instagram',
    # 'allauth.socialaccount.providers.linkedin',
    # 'allauth.socialaccount.providers.linkedin_oauth2',
    # 'allauth.socialaccount.providers.odnoklassniki',
    # 'allauth.socialaccount.providers.openid',
    # 'allauth.socialaccount.providers.persona',
    # 'allauth.socialaccount.providers.soundcloud',
    # 'allauth.socialaccount.providers.stackexchange',
    # 'allauth.socialaccount.providers.tumblr',
    # 'allauth.socialaccount.providers.twitch',
    # 'allauth.socialaccount.providers.twitter',
    # 'allauth.socialaccount.providers.vimeo',
    # 'allauth.socialaccount.providers.vk',
    # 'allauth.socialaccount.providers.weibo',
    # 'allauth.socialaccount.providers.xing',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'biostar3.middleware.GlobalMiddleware',
)

ROOT_URLCONF = 'biostar3.urls'

WSGI_APPLICATION = 'biostar3.wsgi.application'

AUTH_USER_MODEL = 'forum.User'

BIOSTAR_HOME = get_env('BIOSTAR_HOME')

TEMPLATE_PATH = abspath(get_env('THEME_PATH'))
DEFAULT_PATH = abspath(BIOSTAR_HOME, "themes", "default")

TEMPLATE_DIRS = (
    TEMPLATE_PATH,
    DEFAULT_PATH,
)

TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
)

TEMPLATE_STRING_IF_INVALID = "*** missing variable ***"

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
    'biostar3.context.extras',
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
STATIC_ROOT = abspath(BIOSTAR_HOME, "export", "static")

# User uploaded media files
MEDIA_URL = '/media/'
MEDIA_ROOT = abspath(BIOSTAR_HOME, "export", "media")

# The default logo for all groups
DEFAULT_GROUP_LOGO = "images/logo.png"

STATICFILES_DIRS = (
    abspath(TEMPLATE_PATH, "static"),
)

STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    'compressor.finders.CompressorFinder',
)

# The CSS classes associated with the Django messages framework.
MESSAGE_TAGS = {
    10: 'info', 20: 'info', 25: 'success', 30: 'warning', 40: 'error',
}

# Haystack data connection.
HAYSTACK_CONNECTIONS = {
    'default': {
        'ENGINE': 'haystack.backends.whoosh_backend.WhooshEngine',
        'PATH': get_env('SEARCH_INDEX'),
    },
}

AUTHENTICATION_BACKENDS = (
    # Needed to login by username in Django admin, regardless of `allauth`
    "django.contrib.auth.backends.ModelBackend",

    # `allauth` specific authentication methods, such as login by e-mail
    "allauth.account.auth_backends.AuthenticationBackend",
)

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'unique-snowflake'
    }
}


SOCIALACCOUNT_PROVIDERS = {
    'persona': {
        'AUDIENCE': 'http://www.lvh.me:8080/',
        'REQUEST_PARAMETERS': {'siteName': 'Biostars'}
    }
}

