from __future__ import absolute_import
from datetime import timedelta
from celery.schedules import crontab

#CELERY_RESULT_BACKEND = 'djcelery.backends.database:DatabaseBackend'
CELERY_RESULT_BACKEND = 'django-db'

BROKER_URL = 'redis://127.0.0.1:7777'

CELERY_TASK_SERIALIZER = 'pickle'

CELERY_ACCEPT_CONTENT = ['pickle']

#CELERYBEAT_SCHEDULE = {}

CELERY_TIMEZONE = 'UTC'