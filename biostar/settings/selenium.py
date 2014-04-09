# -*- coding: utf8 -*-
#
# Django settings for selenium testing
#
from biostar.settings.base import *

# Turn off captcha
CAPTCHA = False

# How many non top level posts per day for users.
MAX_POSTS_NEW_USER = 100
MAX_POSTS_TRUSTED_USER = 100

# How many top level posts per day for a new user.
MAX_TOP_POSTS_NEW_USER = 100
MAX_TOP_POSTS_TRUSTED_USER = 100

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'unique-snowflake'
    }
}