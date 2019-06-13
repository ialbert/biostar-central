from biostar.engine.settings import *
from biostar.forum.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ENGINE_APPS + ACCOUNTS_APPS + EMAILER_APP

ROOT_URLCONF = 'biostar.test.test_urls'

DEBUG = True

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
