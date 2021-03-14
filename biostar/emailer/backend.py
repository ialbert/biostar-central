import logging
import smtplib

from django.core.mail.backends import smtp
from django.core.mail.utils import DNS_NAME

logger = logging.getLogger("engine")


class SSLEmailBackend(smtp.EmailBackend):
    """
    Required for using Amazon SES"
    """

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
