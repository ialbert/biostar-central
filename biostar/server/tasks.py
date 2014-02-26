from __future__ import absolute_import
from django.conf import settings
from datetime import timedelta
from celery.schedules import crontab

import logging

logger = logging.getLogger(__name__)

from celery import Celery

app = Celery('biostar', broker=settings.BROKER_URL)

@app.task
def add(x, y):
    return x + y

@app.task
def test(*args, **kwds):
    logger.info("*** delayed task %s, %s" % (args, kwds))
    return 1000


CELERYBEAT_SCHEDULE = {
    # Executes every Monday morning at 7:30 A.M
    'add-every-monday-morning': {
        'task': 'tasks.add',
        'schedule': crontab(hour=7, minute=30, day_of_week=1),
        'args': (16, 16),
    },
}