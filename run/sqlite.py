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

# Needs to match the server domain.
SESSION_COOKIE_DOMAIN = ".lvh.me"

# This must be set correctly in production.
ALLOWED_HOSTS = [".lvh.me"]

# The secret key can be used to log into the admin account!
# Make sure to change it in production.
SECRET_KEY = '1@localhost.com'

# The database name must be present.
DATABASE_NAME = get_env('DATABASE_NAME')

# Set up the databases.
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': DATABASE_NAME,
    }
}

