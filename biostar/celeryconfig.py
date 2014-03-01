from __future__ import absolute_import
from datetime import timedelta
from celery.schedules import crontab

CELERY_RESULT_BACKEND = 'djcelery.backends.database:DatabaseBackend'

BROKER_URL = 'django://'

CELERY_TASK_SERIALIZER = 'pickle'

CELERY_ACCEPT_CONTENT = ['pickle' ]

CELERYBEAT_SCHEDULE = {
    'test_task': {
        'task': 'biostar.celery.test',
        'schedule': timedelta(seconds=5),
        'args': (16, 16),
    },
}

CELERY_TIMEZONE = 'UTC'