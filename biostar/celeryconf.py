from __future__ import absolute_import
from datetime import timedelta
from celery.schedules import crontab


RESULT_BACKEND = 'django-db'

BROKER_URL = 'redis://127.0.0.1:6379'

# Mitigate deprecation error:
#     The 'BROKER_URL' setting is deprecated and scheduled for removal in
#     version 6.0.0. Use the broker_url instead

broker_url = BROKER_URL

TASK_SERIALIZER = 'json'

ACCEPT_CONTENT = ['json']

BEAT_SCHEDULER = "django_celery_beat.schedulers:DatabaseScheduler"

BEAT_SCHEDULE = {

    'cleanup': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(days=1),
        'args': ["cleanup"],
    },
    'hourly_feed': {
        'task': 'biostar.celery.call_command',
        'schedule': crontab(minute=10),
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
TIMEZONE = 'UTC'