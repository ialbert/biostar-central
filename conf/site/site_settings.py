import logging
from biostar.forum.settings import *
#from biostar.recipes.settings import *

logger = logging.getLogger("biostar")

# Debugging flag.
DEBUG = True

# Set your known secret key here.
SECRET_KEY = "secretkey"

# Admin users will be created automatically with DEFAULT_ADMIN_PASSWORD.
ADMINS = [
    ("Admin User", "admin@localhost")
]

# Set the default admin password.
DEFAULT_ADMIN_PASSWORD = SECRET_KEY

# Set the site domain.
SITE_DOMAIN = "foo.com"

SITE_ID = 1
HTTP_PORT = ''
PROTOCOL = 'http'

ALLOWED_HOSTS = [SITE_DOMAIN]

DATABASE_NAME = "biostardb"

DATABASES = {

    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': DATABASE_NAME,
        'USER': '',
        'PASSWORD': '',
        'HOST': '/var/run/postgresql/',
        'PORT': '',
    },
}

WSGI_APPLICATION = 'conf.run.site_wsgi.application'

# Valid options; block, disable, threaded, uwsgi, celery.
TASK_RUNNER = 'block'

SESSION_COOKIE_SECURE = True

EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'

try:
    # Attempts to load site secrets.
    from .site_secrets import *

    logger.info("Imported settings from .site_secrets")
except ImportError as exc:
    logger.warn(f"No secrets module could be imported: {exc}")

logger.debug(f"SITE_DOMAIN: {SITE_DOMAIN}")
