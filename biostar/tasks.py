from __future__ import absolute_import
from django.conf import settings
from datetime import timedelta
from celery.schedules import crontab
from django.core.mail import get_connection
from celery.utils.log import get_task_logger

logger = get_task_logger(__name__)

from celery import Celery

app = Celery('biostar',
             broker=settings.BROKER_URL)

@app.task
def test(*args, **kwds):
    logger.info("*** delayed task %s, %s" % (args, kwds))
    return 1000

#
# Based on django-celery-email
# https://github.com/pmclanahan/django-celery-email
#
CONFIG = getattr(settings, 'CELERY_EMAIL_TASK_CONFIG', {})
BACKEND = getattr(settings, 'CELERY_EMAIL_BACKEND',
                  'django.core.mail.backends.smtp.EmailBackend')
TASK_CONFIG = {
    'name': 'tasks.send_email',
    'ignore_result': True,
}
TASK_CONFIG.update(CONFIG)

@app.task(**TASK_CONFIG)
def send_email(message, **kwargs):
    conn = get_connection(backend=BACKEND,
                          **kwargs.pop('_backend_init_kwargs', {}))
    try:
        result = conn.send_messages([message])
        logger.debug("Successfully sent email message to %r.", message.to)
        return result
    except Exception as e:
        print(dir(message))
        logger.error("Error sending email to %r: %s retrying.",
                     message.to, e)
        send_email.retry(exc=e)

CELERYBEAT_SCHEDULE = {
    # Executes every Monday morning at 7:30 A.M
    'add-every-monday-morning': {
        'task': 'tasks.add',
        'schedule': crontab(hour=7, minute=30, day_of_week=1),
        'args': (16, 16),
    },
}
