import os
from django.core.cache import cache
from django.core.servers import basehttp
#
# To see all log messages: export DJANGO_LOG_LEVEL=DEBUG
#
LOG_LEVEL = os.getenv('ENGINE_LOG_LEVEL') or 'INFO'

DJANGO_LOG_LEVEL = os.getenv('DJANGO_LOG_LEVEL') or 'WARNING'


class RateLimitFilter(object):
    """
    Limits the number of error emails when errors get triggered.
    """
    TIMEOUT = 600
    CACHE_KEY = "error-limiter"

    def filter(self, record):

        exists = cache.get(self.CACHE_KEY)
        if not exists:
            cache.set(self.CACHE_KEY, 1, self.TIMEOUT)

        return not exists


LOGGING = {

    'version': 1,

    'disable_existing_loggers': True,

    'formatters': {

        'verbose': {
            'format': '%(levelname)s\t%(asctime)s\t%(module)s.%(funcName)s\t%(lineno)s\t%(message)s\t'
        },

        'simple': {
            'format': '%(levelname)s\t%(module)s\t%(message)s'
        },

    },

    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        },

        'rate_limit': {
            '()': 'biostar.logconf.RateLimitFilter',
        },
    },

    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'verbose'
        },

        'mail_admins': {
            'level': 'ERROR',
            'filters': ['require_debug_false', 'rate_limit'],
            'class': 'django.utils.log.AdminEmailHandler',
            'include_html': True,
        }

    },

    'loggers': {

        'django': {
            'handlers': ['console', 'mail_admins'],
            'level': DJANGO_LOG_LEVEL,

        },

        'engine': {
            'handlers': ['console', 'mail_admins'],
            'level': LOG_LEVEL,
        },

        'biostar': {
            'handlers': ['console', 'mail_admins'],
            'level': LOG_LEVEL,
        },

    },
}
