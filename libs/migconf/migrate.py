from main.settings import *
import os

ADMINS = (
    ('Istvan Albert', 'istvan.albert@gmail.com'),
)

SITE_DOMAIN = 'www.biostars.org'

SECRET_KEY = 'secret-key-goes-here'

# set location relative to the current file directory
__CURR_DIR    = path(os.path.dirname(__file__))
HOME_DIR      = path(__CURR_DIR, '..', 'main' )
TEMPLATE_DIR  = path(HOME_DIR, 'templates')

STATIC_DIR    = path(HOME_DIR, 'static')
EXPORT_DIR    = path(HOME_DIR, '..', 'apache', 'export')
WHOOSH_INDEX  = path(HOME_DIR, 'db', 'index')
STATIC_ROOT   = EXPORT_DIR

GOOGLE_TRACKER = "UA-101522-13"

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'biostar-migrate-database',
        'USER': os.environ['PG_USER'],
        'PASSWORD': os.environ['PG_PASSWD'],
        'HOST': '',
        'PORT': '',
    }
}
