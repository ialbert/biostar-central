# Import all values from the base then override site specific settings.
from biostar3.settings.base import *

# Sqlite specific settings.
DEBUG = True
TEMPLATE_DEBUG = True

# The database name must be present.
DATABASE_NAME = get_env('DATABASE_NAME')

# Set up the databases.
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': DATABASE_NAME,
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
