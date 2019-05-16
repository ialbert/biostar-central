import os

#
# To see all log messages: export DJANGO_LOG_LEVEL=DEBUG
#
LOG_LEVEL = os.getenv('DJANGO_LOG_LEVEL') or 'WARNING'

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

    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'verbose'
        },
    },

    # The valid loggers.
    'loggers': {

        'django': {
            'handlers': ['console'],
            'level': 'INFO',
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
