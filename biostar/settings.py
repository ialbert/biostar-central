import os

# Apply the logger settings.
from biostar.logconf import LOGGING

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

# Django debug flag.
DEBUG = True

# Override compression if needed.
# COMPRESS_ENABLED = True

# Change this in production!
SECRET_KEY = '1234'

# The password for admin users. Must be changed in production.
DEFAULT_ADMIN_PASSWORD = "1234"

DEFAULT_FROM_EMAIL = "admin@localhost"


# Admin users will be created automatically with DEFAULT_ADMIN_PASSWORD.
ADMINS = [
    ("Admin User", "admin@localhost")
]

# Maximum amount of cumulative uploaded files a user is allowed, in mega-bytes.
MAX_UPLOAD_SIZE = 300

# Set these for remote hosts.
SITE_ID = 1
SITE_DOMAIN = "localhost"
SITE_NAME = "Biostar Engine"

# Should the site allow signup.
ALLOW_SIGNUP = False

# Helper function for building absolute paths.
def join(*args):
    return os.path.abspath(os.path.join(*args))

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(join(__file__))

# Application specific parameters.
PROTOCOL = "http"
HTTP_PORT = ':8000'

INSTALLED_APPS = [
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
    'pagedown',

    # The order of apps matters in the template loading.
    'biostar.engine.apps.EngineConfig',
    'biostar.emailer.apps.EmailerConfig',
    'biostar.accounts.apps.AccountsConfig',

    'biostar.ftpserver',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'biostar.engine.middleware.engine_middleware',
]

ROOT_URLCONF = 'biostar.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'APP_DIRS': True,
        'OPTIONS': {
            'string_if_invalid': "**MISSING**",
            'context_processors': [
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.template.context_processors.media',
                'django.contrib.messages.context_processors.messages',
                'biostar.engine.context.engine',

            ],
        },
    },
]

WSGI_APPLICATION = 'biostar.wsgi.application'

# Database
# https://docs.djangoproject.com/en/1.11/ref/settings/#databases

DATABASE_NAME = join(BASE_DIR, '..', 'export', 'database', 'engine.db')

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': DATABASE_NAME,
    }
}

# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },

]

ALLOWED_HOSTS = ['www.lvh.me', 'localhost', '127.0.0.1']

# Internationalization
# https://docs.djangoproject.com/en/1.11/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.11/howto/static-files/

MEDIA_ROOT = join(BASE_DIR, '..', 'export', 'media')

# The location of resusable data.
LOCAL_ROOT = join(BASE_DIR, '..', 'export', 'local')

# The location for the table of contents.
TOC_ROOT = join(MEDIA_ROOT, 'tocs')

# Make the table of contents.
os.makedirs(TOC_ROOT, exist_ok=True)

# Sendfile settings go here.
SENDFILE_ROOT = MEDIA_ROOT
SENDFILE_URL = '/protected/'

#SENDFILE_BACKEND = "sendfile.backends.nginx"
SENDFILE_BACKEND = "sendfile.backends.development"

MEDIA_URL = '/media/'
STATIC_URL = '/static/'
STATIC_ROOT = join(BASE_DIR, '..', 'export', 'static')
STATICFILES_DIRS = [
    join(BASE_DIR, "engine", "static"),
]

STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    'compressor.finders.CompressorFinder',
)

# Apply default logger setting.
LOGGER_NAME = "engine"

# Use django-mailer to store emails in the database.
# EMAIL_BACKEND = "mailer.backend.DbBackend"

# The email delivery engine.
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
