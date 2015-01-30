# -*- coding: utf8 -*-
#
# Setting used when running a Postgres based server.
#
from biostar3.settings.base import *

# Site administrators. Make sure to override this.
ADMINS = (
    ("Biostar Community", "1@localhost.com"),
)

MANAGERS = ADMINS

# Get the secret key from the environment.
SECRET_KEY = get_env("SECRET_KEY")

DATABASE_NAME = get_env('DATABASE_NAME')
DATABASE_USER = get_env('DATABASE_USER')
DATABASE_HOST = get_env('DATABASE_HOST')

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': DATABASE_NAME,
        'USER': DATABASE_USER,
        'PASSWORD': '',
        'HOST':  DATABASE_HOST,
        'PORT': 5432,
     }
}

# What hosts may connect to the site.
ALLOWED_HOSTS = ["localhost"]

# Haystack data connection.
HAYSTACK_CONNECTIONS = {
    'default': {
        'ENGINE': 'haystack.backends.whoosh_backend.WhooshEngine',
        'PATH': get_env('SEARCH_INDEX'),
    },
}
