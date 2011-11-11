#
# default setting for the development server
# modify these settings for any publicly facing website
#
import os, sys, re

# on deployed servers make this unique, and don't share it with anybody.
SECRET_KEY = '007'

# turn off debug mode on deployed servers
DEBUG = True

# template debug mode
TEMPLATE_DEBUG = DEBUG

# admin site may fail if this setting is active
TEMPLATE_STRING_IF_INVALID = "*** MISSING ***"

ADMINS = (
    ('Istvan Albert', 'istvan.albert@gmail.com'),
)

MANAGERS = ADMINS

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

# displays debug comments when the server is run from this IP
INTERNAL_IPS = ('127.0.0.1', )

# the directory that this file is located in
__CURR_DIR = path(os.path.dirname(__file__))

# set location relative to the current file directory
HOME_DIR      = path(__CURR_DIR )
DATABASE_DIR  = path(HOME_DIR, 'db')
DATABASE_NAME = path(DATABASE_DIR, 'biostar.db')
TEMPLATE_DIR  = path(HOME_DIR, 'templates')
STATIC_DIR    = path(HOME_DIR, 'static')
EXPORT_DIR    = path(HOME_DIR, '..', 'export')

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3', # Add 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': DATABASE_NAME,                  # Or path to database file if using sqlite3.
        'USER': '',                      # Not used with sqlite3.
        'PASSWORD': '',                  # Not used with sqlite3.
        'HOST': '',                      # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '',                      # Set to empty string for default. Not used with sqlite3.
    }
}

# add external dependecies
__ZIP_LIBS =  [
    path(__CURR_DIR, '..', 'libs'),
    path(__CURR_DIR, '..', 'libs', 'libraries.zip'),
]
sys.path.extend(__ZIP_LIBS)

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# On Unix systems, a value of None will cause Django to use the same
# timezone as the operating system.
# If running in a Windows environment this must be set to the same as your
# system time zone.
TIME_ZONE = 'America/Chicago'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# added a custom test runner
TEST_RUNNER='biostar.tests.main.BiostarTest'

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale
USE_L10N = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/media/"
MEDIA_ROOT = '/media/'

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://media.lawrence.com/media/", "http://example.com/media/"
MEDIA_URL = EXPORT_DIR

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = EXPORT_DIR

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATIC_URL = '/static/'

# URL prefix for admin static files -- CSS, JavaScript and images.
# Make sure to use a trailing slash.
# Examples: "http://foo.com/static/admin/", "/static/admin/".
ADMIN_MEDIA_PREFIX = '/static/admin/'

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


# ADMIN_PASSWORD_OVERRIDE allows one to log in as any user by using the SECRET_KEY
# this is needed during testing, possbily for troubleshooting problems for deployment
ADMIN_PASSWORD_OVERRIDE = True

# attempt to get the version number from the repository
BIOSTAR_VERSION = os.popen("git log --pretty=format:%h -1").read()
if not re.match(r'^[a-z,0-9]+$', BIOSTAR_VERSION):
    BIOSTAR_VERSION = 'unknown'

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
    'main.middleware.LastVisit',
    'main.middleware.PermissionsMiddleware'
)

ROOT_URLCONF = 'main.urls'

TEMPLATE_DIRS = (
    # Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
    TEMPLATE_DIR,
)

TEMPLATE_CONTEXT_PROCESSORS = (
    "django.contrib.auth.context_processors.auth",
    "django.core.context_processors.debug",
    'django.core.context_processors.request',
    'django.core.context_processors.static',
    'django.contrib.messages.context_processors.messages',
    "main.context.extras",
)

AUTH_PROFILE_MODULE = "server.UserProfile"

ROOT_URLCONF = 'main.urls'

AUTHENTICATION_BACKENDS = (
    'django_openid_auth.auth.OpenIDBackend',
    'django.contrib.auth.backends.ModelBackend',
)

OPENID_CREATE_USERS = True
OPENID_UPDATE_DETAILS_FROM_SREG = True
OPENID_USE_AS_ADMIN_LOGIN = True

LOGIN_URL = '/openid/login/'
LOGIN_REDIRECT_URL = '/'

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.humanize',
    # Uncomment the next line to enable the admin:
    'django.contrib.admin',
    # Uncomment the next line to enable admin documentation:
    'django.contrib.admindocs',
    'main.server',
    'django_openid_auth',
    #'taggit',
)

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'class': 'django.utils.log.AdminEmailHandler'
        }
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
        },
    }
}

# version check, we can do it at the end since
# the version is only required in subsequent modules
if sys.version_info < (2, 6):
    sys.stderr.write( '*** this code requires python 2.6 or higher ***' )
    sys.exit()