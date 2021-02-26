from .settings import *

# Do not multi-thread tests.
MULTI_THREAD = False

# Skip hitting the spam indexe when creating test posts
CLASSIFY_SPAM = False


# Turn the emailing tasks off for tests
SEND_MAIL = False

# Default cache
CACHES = {
    'default': {
        #'BACKEND': 'django.core.cache.backends.dummy.DummyCache',
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'unique-snowflake',
    }
}