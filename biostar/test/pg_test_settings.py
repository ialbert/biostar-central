from biostar.recipes.settings import *
from biostar.forum.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ENGINE_APPS + PLANET_APPS + ACCOUNTS_APPS + EMAILER_APP

ROOT_URLCONF = 'biostar.test.test_urls'

MULTI_THREAD = False
DEBUG = True
INIT_PLANET = False

# reCaptcha left alone during testing
RECAPTCHA_PUBLIC_KEY = ''
RECAPTCHA_PRIVATE_KEY = ''

DATABASE_NAME = os.environ.setdefault("DATABASE_NAME", "database.db")

DATABASES = {

    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': DATABASE_NAME,
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    },
}
