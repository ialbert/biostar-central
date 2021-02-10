from __future__ import absolute_import
from datetime import timedelta
from celery.schedules import crontab

RESULT_BACKEND = 'django-db'

# Mitigate deprecation error:
#     The 'BROKER_URL' setting is deprecated and scheduled for removal in
#     version 6.0.0. Use the broker_url instead

#BROKER_URL = 'redis://127.0.0.1:6379'
broker_url = 'redis://127.0.0.1:6379'

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

    'test': {
        'task': 'biostar.celery.test',
        'schedule': 30,
        'args': ["index"],
        'kwargs': {"update": 1}
    },

    'daily_feed': {
        'task': 'biostar.celery.call_command',
        'schedule': crontab(hour='*/2', minute=15),
        'args': ["planet"],
        'kwargs': {"download": True}
    },

    'update_index': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(minutes=10),
        'args': ["index"],
        'kwargs': {"index": 5000}
    },

    # 'awards': {
    #     'task': 'biostar.celery.call_command',
    #     'schedule': timedelta(hours=3),
    #     'args': ["user_crawl"],
    #     'kwargs': {"award": True}
    # },
    #
    # 'hourly_dump': {
    #     'task': 'biostar.celery.call_command',
    #     'schedule': crontab(minute=10),
    #     'args': ["biostar_pg_dump"],
    #     'kwargs': {"hourly": True}
    # },
    #
    # 'daily_dump': {
    #     'task': 'biostar.celery.call_command',
    #     'schedule': crontab(hour=22),
    #     'args': ["biostar_pg_dump"],
    # },


}
TIMEZONE = 'UTC'