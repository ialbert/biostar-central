from biostar.settings import *

DEBUG = True

WSGI_APPLICATION = 'conf.postgres.postgres_wsgi.application'

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'biostar_engine',
        'USER': 'ialbert',
        'PASSWORD': 'ialbert',
        'HOST': 'localhost',
        'PORT': '',
    }
}

try:
    from .postgres_secrets import *
except ImportError as exc:
    print("No postgres_secrets module could be imported")
