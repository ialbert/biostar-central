from __future__ import absolute_import
from datetime import timedelta
from celery.schedules import crontab

CELERY_RESULT_BACKEND = 'djcelery.backends.database:DatabaseBackend'

BROKER_URL = 'django://'

CELERY_TASK_SERIALIZER = 'pickle'

CELERY_ACCEPT_CONTENT = ['pickle']

CELERYBEAT_SCHEDULE = {

    'prune_data': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(days=1),
        'kwargs': dict(name="prune_data")
    },

    'sitemap': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(hours=6),
        'kwargs': dict(name="sitemap")
    },

    'update_index': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(minutes=15),
        'args': ["update_index"],
        'kwargs': {"age": 1}
    },

    'hourly_dump': {
        'task': 'biostar.celery.call_command',
        'schedule': crontab(minute=10),
        'args': ["biostar_pg_dump"],
        'kwargs': {"hourly": True}
    },

    'daily_dump': {
        'task': 'biostar.celery.call_command',
        'schedule': crontab(hour=22),
        'args': ["biostar_pg_dump"],
    },

}

CELERY_TIMEZONE = 'UTC'