from __future__ import absolute_import
from datetime import timedelta
from celery.schedules import crontab

#CELERY_RESULT_BACKEND = 'djcelery.backends.database:DatabaseBackend'
CELERY_RESULT_BACKEND = 'django-db'

BROKER_URL = 'redis://127.0.0.1:6379'

CELERY_TASK_SERIALIZER = 'pickle'

CELERY_ACCEPT_CONTENT = ['pickle']

CELERYBEAT_SCHEDULER = "django_celery_beat.schedulers:DatabaseScheduler"

CELERYBEAT_SCHEDULE = {

    'cleanup': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(days=1),
        'args': ["cleanup"],
    },
    'hourly_feed': {
        'task': 'biostar.celery.call_command',
        'schedule': crontab(minute='*/2'),
        'args': ["planet"],
        'kwargs': {"update": 1}
    },

    'daily_feed': {
        'task': 'biostar.celery.call_command',
        'schedule': crontab(hour='*/2', minute=15),
        'args': ["planet"],
        'kwargs': {"download": True}
    },

}
CELERY_TIMEZONE = 'UTC'