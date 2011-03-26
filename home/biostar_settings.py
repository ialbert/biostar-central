# Django settings for biostar project.
import os, sys, random, hashlib

# turn this off in production
DEBUG = True
TEMPLATE_DEBUG = DEBUG

def path_join(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

# the directory that this file is located in
__CURR_DIR = path_join(os.path.dirname(__file__))

# some dependecies may be distributed as zipfiles
__ZIP_LIBS =  [
    path_join(__CURR_DIR, '..', 'biostar', 'libs', 'openid-libraries.zip'),
]
sys.path.extend(__ZIP_LIBS)

# set paths to various file locations
HOME_DIR = path_join(__CURR_DIR )
DATABASE_DIR = path_join(HOME_DIR, 'db')
DATABASE_NAME = path_join(DATABASE_DIR, 'biostar.db')
TEMPLATE_DIR = path_join(HOME_DIR, 'templates')
STATIC_DIR = path_join(HOME_DIR, 'static')

ADMINS = (
    ('Istvan Albert', 'istvan.albert@gmail.com'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3', # Add 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': DATABASE_NAME,                      # Or path to database file if using sqlite3.
        'USER': '',                      # Not used with sqlite3.
        'PASSWORD': '',                  # Not used with sqlite3.
        'HOST': '',                      # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '',                      # Set to empty string for default. Not used with sqlite3.
    }
}

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

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale
USE_L10N = True

# Absolute path to the directory that holds media.
# Example: "/home/media/media.lawrence.com/"
MEDIA_ROOT = ''

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash if there is a path component (optional in other cases).
# Examples: "http://media.lawrence.com", "http://example.com/media/"
MEDIA_URL = ''

# URL prefix for admin media -- CSS, JavaScript and images. Make sure to use a
# trailing slash.
# Examples: "http://foo.com/media/", "/media/".
ADMIN_MEDIA_PREFIX = '/media/'

# On first run generates a unique secret key that also serves as the administrative password
SECRET_FILE = path_join(HOME_DIR, 'secret-key.txt')
if not os.path.isfile(SECRET_FILE):
    fp = file(SECRET_FILE, 'wt')
    val = str(random.getrandbits(128))
    val = hashlib.md5(val).hexdigest()
    fp.write(val)
    fp.close()

SECRET_KEY = file(SECRET_FILE).read().strip()

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

TEMPLATE_STRING_IF_INVALID = "*** MISSING ***"

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    
)

AUTH_PROFILE_MODULE = "server.UserProfile"

ROOT_URLCONF = 'biostar.urls'

TEMPLATE_DIRS = (
    TEMPLATE_DIR,
)

AUTHENTICATION_BACKENDS = (
    'django_openid_auth.auth.OpenIDBackend',
    'django.contrib.auth.backends.ModelBackend',
)

OPENID_CREATE_USERS = True
OPENID_UPDATE_DETAILS_FROM_SREG = True

LOGIN_URL = '/openid/login/'
LOGIN_REDIRECT_URL = '/'

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    # Uncomment the next line to enable the admin:
    'django.contrib.admin',
    # Uncomment the next line to enable admin documentation:
    # 'django.contrib.admindocs',
    'biostar.server',
    'django_openid_auth',
)

# version check
if sys.version_info < (2, 6):
    sys.stderr.write( '*** this code requires python 2.6 or higher ***' )
    sys.exit()
