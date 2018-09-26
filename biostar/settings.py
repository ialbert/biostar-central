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
SECRET_KEY = 'secret-key'

# Private key used to validate external logins
LOGIN_PRIVATE_KEY = "private-key"

# The password for admin users. Must be changed in production.
DEFAULT_ADMIN_PASSWORD = "admin@localhost"

RECAPTCHA_PUBLIC_KEY = ""
RECAPTCHA_PRIVATE_KEY = ""

POSTS_PER_PAGE = 35
USERS_PER_PAGE = 100
MESSAGES_PER_PAGE = 100
TAGS_PER_PAGE = 50

VOTE_FEED_COUNT = 10
LOCATION_FEED_COUNT = 5
AWARDS_FEED_COUNT = 10
REPLIES_FEED_COUNT = 15

ENGINE_AS_ROOT = True

SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True

# Set the home page to the engine or forum
INTERNAL_IPS = ['127.0.0.1']

# Admin users will be created automatically with DEFAULT_ADMIN_PASSWORD.
ADMINS = [
    ("Admin User", "admin@localhost")
]

# The default sender name on emails.
DEFAULT_FROM_EMAIL = f"{ADMINS[0][0]} <{ADMINS[0][1]}>"

# Maximum amount of cumulative uploaded files a user is allowed, in mega-bytes.
MAX_UPLOAD_SIZE = 10

# These must be set remote hosts.
SITE_ID = 1
SITE_DOMAIN = "localhost"
SITE_NAME = "Biostar Engine"

# Deployment specific parameters.
PROTOCOL = "http"
HTTP_PORT = ':8000'

FTP_HOST = "localhost"
FTP_PORT = 8021

# Should the site allow signup.
ALLOW_SIGNUP = True

# Allow users to toggle their status moderator
ALLOW_SELF_MODERATE = False

# Maximum size of each file upload in MB
MAX_FILE_SIZE_MB = 300

LOGIN_REDIRECT_URL = "/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

# Helper function for building absolute paths.
def join(*args):
    return os.path.abspath(os.path.join(*args))


# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(join(__file__))

SOCIALACCOUNT_ADAPTER = "biostar.accounts.adapter.SocialAccountAdapter"

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
    'taggit',
    'debug_toolbar',
    'snowpenguin.django.recaptcha2',

    # The order of apps matters in the template loading
    'biostar.engine.apps.EngineConfig',
    'biostar.emailer.apps.EmailerConfig',
    'biostar.accounts.apps.AccountsConfig',
    'biostar.forum.apps.ForumConfig',
    'biostar.message.apps.MessageConfig',
    'biostar.ftpserver',

    # Allauth templates come last.
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
    'allauth.socialaccount.providers.github',
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
    'debug_toolbar.middleware.DebugToolbarMiddleware'

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

AUTHENTICATION_BACKENDS = (
    "django.contrib.auth.backends.ModelBackend",
    "allauth.account.auth_backends.AuthenticationBackend",
)

WSGI_APPLICATION = 'biostar.wsgi.application'

# Database settings.
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

# The location of application specific data.
LOCAL_ROOT = join(BASE_DIR, '..', 'export', 'local')

# The location for the table of contents.
TOC_ROOT = join(MEDIA_ROOT, 'tocs')

# Time between two accesses from the same IP to qualify as a different view.
POST_VIEW_MINUTES = 7

# Configure language detection
LANGUAGE_DETECTION = ['en']

# Ensure that the table of directory exists.
os.makedirs(TOC_ROOT, exist_ok=True)

# Sendfile settings go here.
SENDFILE_ROOT = MEDIA_ROOT
SENDFILE_URL = '/protected/'

# Settings used to enable/disable the forum.
ENABLE_FORUM = True
ONLY_FORUM_URLS = False

# Session settings go here
SESSION_KEY = "session"

COUNT_INTERVAL_WEEKS = 10000

SESSION_ENGINE = "django.contrib.sessions.backends.cached_db"

# SENDFILE_BACKEND = "sendfile.backends.nginx"
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

# List of social login clients tuples.
# ( name, client_id, secret )

LOGIN_REDIRECT_URL = "/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

# Default clients redirect to localhost.
# Default clients may not be operational. See the
# django allauth documentataion on how to set them up.

#
# Callback example settings:
#
# http://localhost:8000/accounts/social/google/login/callback/
# http://localhost:8000/accounts/social/github/login/callback/
#
SOCIAL_CLIENTS = [

    ("Google", "547073349197-ri3eku9fdpi1ble7eoc4amlsrh7m2oiv.apps.googleusercontent.com", "DR1-zMqOLTqRGvhSDb5rQBMg"),
    ("GitHub", "d8493ce8967ea5abbd73", "04ff5043ebce2317d3c26bfe90fbc5e67fa38d05")

]
