from __future__ import absolute_import
from django.conf import settings
import functools
from celery.utils.log import get_task_logger
from celery.schedules import crontab
import os

#from biostar import const
from datetime import timedelta

logger = get_task_logger(__name__)

from celery import Celery

# Django settings module
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'biostar.settings')

app = Celery('biostar')

# Read the configuration from the config file.
app.config_from_object('biostar.celeryconf')

@app.task
def call_command(name, *args, **kwargs):
    "Calls a django command in a delayed fashion"
    logger.info(f"calling django command {name} with {args} and {kwargs}")
    from django.core.management import call_command
    call_command(name, *args, **kwargs)

@app.task
def test(*args, **kwds):
    logger.info(f"*** executing task {__name__} {args}, {kwds}")

