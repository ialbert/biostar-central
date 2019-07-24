import os

#
# To see all log messages: export DJANGO_LOG_LEVEL=DEBUG
#
LOG_LEVEL = os.getenv('DJANGO_LOG_LEVEL') or 'INFO'

DJANGO_LOG = 'WARNING'


class RateLimitFilter(object):
    """
    Limits the number of error emails when errors get triggered.
    """

    def filter(self, record):
        from django.core.cache import cache
        TIMEOUT = 600
        CACHE_KEY = "error-limiter"

        exists = cache.get(CACHE_KEY)
        if not exists:
            cache.set(CACHE_KEY, 1, TIMEOUT)

        return not exists


LOGGING = {

    'version': 1,

    'disable_existing_loggers': False,

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
            'handlers': ['console'],
            'level': DJANGO_LOG,
        },

        'engine': {
            'handlers': ['console'],
            'level': LOG_LEVEL,
        },

        'biostar': {
            'handlers': ['console'],
            'level': LOG_LEVEL,
        },

    },
}
