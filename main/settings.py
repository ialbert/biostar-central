#
# default setting for the development server
# modify these settings for any publicly facing website
#
import os, sys, re

# database migrations via Django South
import south

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))


# on deployed servers make this unique, and don't share it with anybody.
SECRET_KEY = '007'

# turn off debug mode on deployed servers
DEBUG = True

# template debug mode
TEMPLATE_DEBUG = DEBUG

# admin site may fail if this setting is active
TEMPLATE_STRING_IF_INVALID = "*** MISSING ***"

# the time in seconds between session updates
SESSION_UPDATE_TIME = 10 * 60  # in seconds

ADMINS = (
    ('Default Admin', 'your-mail-here@your-server-here.com'),
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
TEMPLATE_DIR  = path(HOME_DIR, 'main', 'templates')
STATIC_DIR    = path(HOME_DIR, 'static')
EXPORT_DIR    = path(HOME_DIR, '..', 'apache', 'export')
WHOOSH_INDEX  = path(HOME_DIR, 'db', 'index')
PLANET_DIR    = path(HOME_DIR, 'db', 'planet')

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

# email specific settings
EMAIL_HOST = 'smtp.yourserver.com'
EMAIL_HOST_USER = 'user'
EMAIL_HOST_PASSWORD = 'password'
DEFAULT_FROM_EMAIL = 'default'
SERVER_EMAIL = 'default'

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

# the dommain for this site
SITE_DOMAIN = 'localhost:8080'

# added a custom test runner
TEST_RUNNER='server.tests.runner.BiostarTest'

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale
USE_L10N = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/media/"
MEDIA_ROOT = "."

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://media.lawrence.com/media/", "http://example.com/media/"
MEDIA_URL = "/"

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = EXPORT_DIR

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
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

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    #'django.middleware.locale.LocaleMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'main.middleware.LastVisit',
)

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.dummy.DummyCache' if DEBUG else 'django.core.cache.backends.locmem.LocMemCache',
        #'BACKEND':  'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'unique-snowflake'
    }
}

#SESSION_ENGINE = "django.contrib.sessions.backends.cache" 

ROOT_URLCONF = 'main.urls'

TEMPLATE_DIRS = [
    # Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
    TEMPLATE_DIR,
]

TEMPLATE_CONTEXT_PROCESSORS = (
    "django.contrib.auth.context_processors.auth",
    "django.core.context_processors.debug",
    'django.core.context_processors.request',
    'django.core.context_processors.static',
    'django.contrib.messages.context_processors.messages',
    'django.core.context_processors.i18n',
    "main.context.extras",
    "main.context.popular_tags"
)

AUTH_PROFILE_MODULE = "server.UserProfile"

ROOT_URLCONF = 'main.urls'

AUTHENTICATION_BACKENDS = (
    'django_openid_auth.auth.OpenIDBackend',
    'django.contrib.auth.backends.ModelBackend',
)

OPENID_CREATE_USERS = True
OPENID_UPDATE_DETAILS_FROM_SREG = False
OPENID_USE_AS_ADMIN_LOGIN = True

# allow migration based on user email
ALLOW_OPENID_MIGRATION = True

LOGIN_URL = '/openid/login/'
LOGIN_REDIRECT_URL = '/'

INSTALLED_APPS = [
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.humanize',
    'django.contrib.markup',
    'django.contrib.messages',
    # Uncomment the next line to enable the admin:
    'django.contrib.admin',
    # Uncomment the next line to enable admin documentation:
    'django.contrib.admindocs',
    'south',
    'main.server',
    'django_openid_auth',
    'django.contrib.sitemaps',
]

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'simple': {
            'format': '%(levelname)s\t%(name)s\t%(funcName)s\t%(asctime)s \t%(message)s'
        },
    },
    'handlers': {
        'console':{
            'level':'DEBUG',
            'class':'logging.StreamHandler',
            'formatter': 'simple'
        },        
    },
    'loggers': {
        'main.server.views': {
            'handlers': [ 'console' ],
            'level': 'INFO',
        },
        'main.server.action': {
            'handlers': [ 'console' ],
            'level': 'INFO',
        },
        'main.server.auth': {
            'handlers': [ 'console' ],
            'level': 'INFO',
        },
        'main.server.models': {
            'handlers': [ 'console' ],
            'level': 'INFO',
        },
        'main.server.search': {
            'handlers': [ 'console' ],
            'level': 'INFO',
        }
    }
}

# google analytics tracker and domain
GOOGLE_TRACKER = ""
GOOGLE_DOMAIN  = ""

# needs to be turned on explicitly
CONTENT_INDEXING = True

# rank gains expressed in hours
POST_UPVOTE_RANK_GAIN = 1
POST_VIEW_RANK_GAIN = 0.1
BLOG_VIEW_RANK_GAIN = 0.1

# if this is set together with the DEBUG mode allows test logins
# don't turn it on in production servers!
SELENIUM_TEST_LOGIN_TOKEN = None

# setting the session for multiple servers
SESSION_COOKIE_DOMAIN = ""

#
# TEMPLATE LAYOUT,
# One may override these variables from the settings file
# 

# this data governs the layout of the PILL_BAR    
# bar name, link url, link name, counter key
USER_PILL_BAR = [
    ("all", "/", "Show&nbsp;All", "" ),
    ("mytags", "/show/mytags/", "My&nbsp;Tags", "" ),
    ("news", "/show/news/", "News", "News" ),
    ("questions", "/show/questions/", "Questions", "Question" ),
    ("unanswered", "/show/unanswered/", "Unanswered", "Unanswered" ),
    ("tutorials", "/show/tutorials/", "Tutorials", "Tutorial" ),
    ("tools", "/show/tools/", "Tools", "Tool" ),
    ("videos", "/show/videos/", "Videos", "Video" ),
    ("jobs", "/show/jobs/", "Jobs", "Job" ),
]

ANON_PILL_BAR = [
    ("all", "/", "Show&nbsp;All", "" ),
    ("news", "/show/news/", "News", "News" ),
    ("questions", "/show/questions/", "Questions", "Question" ),
    ("unanswered", "/show/unanswered/", "Unanswered", "Unanswered" ),
    ("tutorials", "/show/tutorials/", "Tutorials", "Tutorial" ),
    ("tools", "/show/tools/", "Tools", "Tool" ),
    ("videos", "/show/videos/", "Videos", "Video" ),
    ("jobs", "/show/jobs/", "Jobs", "Job" ),
]

#
# remapping the templates to local versions
# a row is the way a post is rendered on a page
# list below the templates to be loaded for a post type
# to reduce clutter there is a default mapper that
# for missing types attempts to map each type to rows/row.type.html
# django template lookup rules apply
#
TEMPLATE_ROWS = {
    'job': "rows/row.job.html",
}

# version check, we can do it at the end since
# the version is only required in subsequent modules
if sys.version_info < (2, 6):
    sys.stderr.write( '*** this code requires python 2.6 or higher ***' )
    sys.exit()