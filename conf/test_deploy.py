#
# import from the main settings then override some of the values
#
from main.settings import *

# set the domain administrators
ADMINS = (
    ('Istvan Albert', 'istvan.albert@gmail.com'),
)

DEBUG = True
# template debug mode
TEMPLATE_DEBUG = DEBUG

# set the site url
SITE_DOMAIN = 'localhost:8080'

# set the cookie domain as needed
SESSION_COOKIE_DOMAIN = ".biostars.org"

# set the secret key for the site
SECRET_KEY = 'murjkj468712u7u2888271209239929u7u2888271209239929u7u28882'

SESSION_UPDATE_TIME = 30

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'biostar-test-database',
        'USER': 'ialbert',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
        }
}

CACHES = {
    'default': {
        'BACKEND':  'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'unique-snowflake'
    }
}

EMAIL_HOST = 'smtp.psu.edu'
EMAIL_HOST_USER = ''
EMAIL_HOST_PASSWORD = ''
DEFAULT_FROM_EMAIL = 'istvan.albert@gmail.com'
SERVER_EMAIL = 'admin@biostars.org'
