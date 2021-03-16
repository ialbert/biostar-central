from __future__ import absolute_import
from datetime import timedelta
from celery.schedules import crontab
from .celery import app
from django.conf import settings

RESULT_BACKEND = 'django-db'

# Mitigate deprecation error:
#     The 'BROKER_URL' setting is deprecated and scheduled for removal in
#     version 6.0.0. Use the broker_url instead

#BROKER_URL = 'redis://127.0.0.1:6379'
broker_url = 'redis://127.0.0.1:6379'

TASK_SERIALIZER = 'json'

ACCEPT_CONTENT = ['json']

BEAT_SCHEDULER = "django_celery_beat.schedulers:DatabaseScheduler"

#app.conf.timezone = 'UTC'
TIMEZONE = 'UTC'

app.conf.broker_url = 'redis://127.0.0.1:6379'

# Discover tasks in applications.
app.autodiscover_tasks(
    lambda: settings.TASK_MODULES
)
