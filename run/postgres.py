# -*- coding: utf8 -*-
#
# Setting used when running a Postgres based server.
#
from biostar3.settings.base import *

# Turn this off during deployment.
DEBUG = True
TEMPLATE_DEBUG = True

# Site administrators. Make sure to override this.
ADMINS = (
    ("Biostar Community", "1@localhost.com"),
)

DEFAULT_FROM_EMAIL = "Site Admin <1@localhost.com>"

MANAGERS = ADMINS

# The secret key can be used to log into the admin email!
# Make sure to change it in production.
SECRET_KEY = 'secret_key'

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

# How many posts per page.
POSTS_PER_PAGE = 35