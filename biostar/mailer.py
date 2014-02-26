import smtplib
import logging

from django.conf import settings
from django.core.mail.utils import DNS_NAME
from django.core.mail.backends import smtp
from django.core.mail.backends.base import BaseEmailBackend

logger = logging.getLogger(__name__)

from biostar.tasks import send_email

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
        results = []
        kwargs['_backend_init_kwargs'] = self.init_kwargs
        for msg in email_messages:
            results.append(send_email.delay(msg, **kwargs))
        return results