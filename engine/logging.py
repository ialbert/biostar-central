LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,

    # What to output for a logging message.
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s %(process)d %(thread)d %(message)s'
        },
        'simple': {
            'format': '%(levelname)s %(module)s.%(funcName)s %(message)s'
        },
    },

    # Destination of the logging messages.
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        },
    },

    # The logger names.
    'loggers': {
        'django': {
            'handlers': ['console'],
            'level': 'ERROR',
        },
        'engine': {
            'handlers': ['console'],
            'level': 'DEBUG',
        },
    },
}
