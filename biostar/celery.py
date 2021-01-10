from __future__ import absolute_import
from django.conf import settings
from celery.utils.log import get_task_logger
import os

#from biostar import const
from datetime import timedelta

logger = get_task_logger(__name__)

from celery import Celery

# Django settings module
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'conf.run.site_settings')

app = Celery('biostar')

# Read the configuration from the config file.
app.config_from_object('biostar.celeryconf')

# Discover tasks in applications.
app.autodiscover_tasks(
    lambda: ["biostar.forum.tasks"]
)

@app.task
def test(*args, **kwds):
    logger.info(f"*** executing task {__name__} {args}, {kwds}")
