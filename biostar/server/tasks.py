from __future__ import absolute_import
from django.conf import settings
from datetime import timedelta
from celery.schedules import crontab
from django.core.mail.backends.base import BaseEmailBackend
from django.core.mail import get_connection
from celery.utils.log import get_task_logger

logger = get_task_logger(__name__)

from celery import Celery

app = Celery('biostar', broker=settings.BROKER_URL)

@app.task
def add(x, y):
    return x + y

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
        logger.warning("Failed to send email message to %r, retrying.",
                       message.to)
        send_email.retry(exc=e)

class CeleryEmailBackend(BaseEmailBackend):
    def __init__(self, fail_silently=False, **kwargs):
        super(CeleryEmailBackend, self).__init__(fail_silently)
        self.init_kwargs = kwargs

    def send_messages(self, email_messages, **kwargs):
        results = []
        kwargs['_backend_init_kwargs'] = self.init_kwargs
        for msg in email_messages:
            results.append(send_email.delay(msg, **kwargs))
        return results

CELERYBEAT_SCHEDULE = {
    # Executes every Monday morning at 7:30 A.M
    'add-every-monday-morning': {
        'task': 'tasks.add',
        'schedule': crontab(hour=7, minute=30, day_of_week=1),
        'args': (16, 16),
    },
}