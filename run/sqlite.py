# Import all values from the base then override site specific settings.
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
