# -*- coding: utf-8 -*-
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

# The secret key can be used to log into the admin account!
# Make sure to change it in production.
SECRET_KEY = '1@localhost.com'

# Needs to match the server domain.
SESSION_COOKIE_DOMAIN = ".lvh.me"

# This must be set correctly in production.
ALLOWED_HOSTS = [".lvh.me"]

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

# How many posts per page.
POSTS_PER_PAGE = 35

# Google Analytics Property ID.
# Please remove the "False" and add your key in between "".
GOOGLE_ANALYTICS_PROPERTY_ID = False