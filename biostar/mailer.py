from __future__ import absolute_import
import smtplib
import logging
from django.conf import settings
from django.core.mail.utils import DNS_NAME
from django.core.mail.backends import smtp
from django.core.mail.backends.base import BaseEmailBackend
from django.core.mail import get_connection

from .celery import app

#logger = logging.getLogger(__name__)

from celery.utils.log import get_task_logger
logger = get_task_logger(__name__)

# Based on django-celery-email
# https://github.com/pmclanahan/django-celery-email
CONFIG = getattr(settings, 'CELERY_EMAIL_TASK_CONFIG', {})

# Get the email sending backend.
BACKEND = getattr(settings, 'CELERY_EMAIL_BACKEND',
                  'django.core.mail.backends.smtp.EmailBackend')

TASK_CONFIG = {
    'name': 'celery.send_email',
    'ignore_result': True,
}
TASK_CONFIG.update(CONFIG)

@app.task
def send_email(message, **kwargs):

    conn = get_connection(backend=BACKEND,
                          **kwargs.pop('_backend_init_kwargs', {}))
    try:
        logger.info("send_email to %s", message.to)
        result = conn.send_messages([message])
        return result
    except Exception as e:
        logger.error("Error sending email from %s to %r: %s",
                     message.from_email, message.to, e)
        #send_email.retry(exc=e)

class SSLEmailBackend(smtp.EmailBackend):
    "Required for Amazon SES"
    def __init__(self, *args, **kwargs):
      kwargs.setdefault('timeout', 5)
      super(SSLEmailBackend, self).__init__(*args, **kwargs)

    def open(self):
        if self.connection:
            return False
        try:
            logger.info("sending email via %s" % self.host)
            self.connection = smtplib.SMTP_SSL(self.host, self.port,
                                               local_hostname=DNS_NAME.get_fqdn())
            if self.username and self.password:
                self.connection.login(self.username, self.password)
            return True
        except:
            if not self.fail_silently:
                raise

class CeleryEmailBackend(BaseEmailBackend):
    def __init__(self, fail_silently=False, **kwargs):
        super(CeleryEmailBackend, self).__init__(fail_silently)
        self.init_kwargs = kwargs

    def send_messages(self, email_messages, **kwargs):
        logger.debug("send_messages %s" % len(email_messages))
        results = []
        kwargs['_backend_init_kwargs'] = self.init_kwargs
        for msg in email_messages:
            results.append(send_email.delay(msg, **kwargs))
        return results