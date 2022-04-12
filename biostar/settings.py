"""
Common settings that are applicable to all biostar apps.
"""

import os

# The logging configuration
from biostar.logconf import LOGGING


# Helper function for building absolute paths.
def join(*args):
    return os.path.abspath(os.path.join(*args))


# Pagedown
PAGEDOWN_IMAGE_UPLOAD_ENABLED = False

LANGUAGE_DETECTION = ["en"]

# Set the home page to the engine or forum
INTERNAL_IPS = ['127.0.0.1']

# Admin users will be created automatically with DEFAULT_ADMIN_PASSWORD.
ADMINS = [
    ("Admin User", "admin@localhost"),
]

DEFAULT_ADMIN_PASSWORD = "admin@localhost"

# Allowed CORS websites
CORS_ORIGIN_WHITELIST = []

POSTGRES_HOST = os.environ.setdefault("POSTGRES_HOST", "")

# Shortcut to first admin information.
ADMIN_NAME, ADMIN_EMAIL = ADMINS[0]

# The default sender name on emails.
DEFAULT_FROM_EMAIL = ADMIN_EMAIL

# The default no reply email.
DEFAULT_NOREPLY_EMAIL = ADMIN_EMAIL

# Email used to send errors to ADMINS
SERVER_EMAIL = ADMIN_EMAIL

FROM_EMAIL_PATTERN = "%s <%s>"

SEND_MAIL = True

# Show debug toolbar
DEBUG_TOOLBAR = False

# The current directory path.
__CURR_DIR = os.path.dirname(join(__file__))

# The directory relative to which all content is stored.
BASE_DIR = join(__CURR_DIR, "..")

# Django debug flag.
DEBUG = True

# Default installed apps.
DEFAULT_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.humanize',
    'mailer',
    'compressor',
    'taggit',
    'snowpenguin.django.recaptcha2',

]

# Enabled apps.
INSTALLED_APPS = DEFAULT_APPS

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',

]

# Site specific information
SITE_ID = 1
SITE_DOMAIN = "localhost"
SITE_NAME = "Biostar Engine"

# Deployment specific parameters.
PROTOCOL = "http"
HTTP_PORT = '8000'
BASE_URL = f"{PROTOCOL}://{SITE_DOMAIN}:{HTTP_PORT}"

# Change this in production!
SECRET_KEY = 'secret-key'

# Change this in production!
API_KEY = "api-key"

# Used during testing
#TEMPLATE_DIRS = [os.path.join(BASE_DIR, 'themes', 'biostar')]

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'APP_DIRS': True,
        #'DIRS': TEMPLATE_DIRS,
        'OPTIONS': {
            'string_if_invalid': "**MISSING**",
            'context_processors': [
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.template.context_processors.media',
                'django.contrib.messages.context_processors.messages',
                'biostar.context.main',
            ],
        },

    },
]

# Authentication backend.
AUTHENTICATION_BACKENDS = [
    "django.contrib.auth.backends.ModelBackend",
]

# The WSGI application.
WSGI_APPLICATION = 'biostar.wsgi.application'

# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators
AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },

]

# Database directory.
# https://docs.djangoproject.com/en/1.11/ref/settings/#databases
DATABASE_DIR = os.path.join(BASE_DIR, 'export', 'db')

DEFAULT_AUTO_FIELD = 'django.db.models.AutoField'

os.makedirs(DATABASE_DIR, exist_ok=True)

DATABASE_NAME = os.environ.setdefault("DATABASE_NAME", "database.db")
# Ensure database is inside database directory.
DATABASE_NAME = os.path.join(DATABASE_DIR, DATABASE_NAME)

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': DATABASE_NAME,
    }
}

ALLOWED_HOSTS = ['www.lvh.me', 'localhost', '127.0.0.1']

# The URL configuration.
ROOT_URLCONF = 'biostar.urls'

# Internationalization
# https://docs.djangoproject.com/en/1.11/topics/i18n/
LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'America/New_York'

USE_I18N = True

USE_L10N = True

USE_TZ = True

# The static URL start.
STATIC_URL = '/static/'

# The static root directory.
STATIC_ROOT = join(BASE_DIR, 'export', 'static')

# Global directories for static files.
STATICFILES_DIRS = [
    join(BASE_DIR, "biostar", "static"),
]

# The media URL start.
MEDIA_URL = '/media/'

# The media root directory.
MEDIA_ROOT = join(BASE_DIR, 'export', 'media')

# The root for all docs
DOCS_ROOT = join(BASE_DIR, 'docs')

STATICFILES_FINDERS = [
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    'compressor.finders.CompressorFinder'
]

# Apply default logger setting.
LOGGER_NAME = "biostar"

# Valid options; block, disabled, threaded, uwsgi, celery.
TASK_RUNNER = 'block'

TASK_MODULES = []

# The email delivery engine.
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'

# Session engine.
SESSION_ENGINE = "django.contrib.sessions.backends.db"

# Session key name.
SESSION_KEY = "session"

# Session key to keep track of counts
SESSION_COUNT_KEY = "counts"

DROPDOWN_TAGS = False
