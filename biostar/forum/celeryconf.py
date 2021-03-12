from datetime import timedelta
from celery.schedules import crontab


BEAT_TASKS = {
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
    'update_index': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(minutes=10),
        'args': ["index"],
        'kwargs': {"index": 5000}
    },

    'awards': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(hours=1),
        'args': ["tasks"],
        'kwargs': {"action": 'award'}
    },

    'hourly_dump': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(hours=1),
        'args': ["tasks"],
        'kwargs': {"action": 'pg_dump', 'hourly': True}
    },
    'daily_dump': {
        'task': 'biostar.celery.call_command',
        'schedule': timedelta(hours=22),
        'args': ["tasks"],
        'kwargs': {"action": 'pg_dump'}
    },
}


TASK_MODULES = ("biostar.forum.tasks", )
