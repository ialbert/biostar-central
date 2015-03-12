"""
A simplifed version of https://github.com/BradWhittington/django-templated-email

Allows an email to be specified in a template.

Extracts all information from a single template like so. All three blocks must be present.

    {% block subject %} declares the subject
    {% block text %} declares text/plain
    {% block html %} declares text/html

"""
import logging, smtplib, string

from django.conf import settings
from django.template.loader import get_template
from django.template.loader_tags import BlockNode
from django.core.mail import EmailMultiAlternatives
from django.template import Context
from django.core.mail.utils import DNS_NAME
from django.core.mail.backends import smtp
from django.core.mail.backends.base import BaseEmailBackend
from django.core.mail import get_connection

logger = logging.getLogger("biostar")


def get_node(template, name):
    for node in template:
        if isinstance(node, BlockNode) and node.name == name:
            return node
    raise Exception("Node '%s' could not be found in template." % name)


class EmailTemplate(object):
    """

    """
    SUBJECT, TEXT, HTML = "subject", "text", "html"

    def __init__(self, template_name, data={}, **kwargs):
        data.update(kwargs)
        context = Context(data, autoescape=False)
        template = get_template(template_name)

        # Extract the blocks for each part of the email.
        subj = get_node(template, name=self.SUBJECT).render(context)
        # Email subject may not contain newlines
        lines = subj.splitlines()
        lines = map(string.strip, lines)

        self.subj = "".join(lines)
        self.text = get_node(template, name=self.TEXT).render(context)
        self.html = get_node(template, name=self.HTML).render(context)


    def send(self, to, from_email=None, cc=None, bcc=None, headers={}, token=None):

        # Support address in the `to` parameter.
        if type(to) != list:
            to = [to]

        if token:
            # Add reply to token to header.
            reply_to = settings.REPLY_TO_PATTERN % token
            headers.update({'Reply-To': reply_to})

        # Fall back to defaults
        from_email = from_email or settings.DEFAULT_FROM_EMAIL
        msg = EmailMultiAlternatives(
            subject=self.subj, body=self.text, from_email=from_email, to=to,
            cc=cc, bcc=bcc, headers=headers,
        )
        if self.html:
            msg.attach_alternative(self.html, 'text/html')

        msg.send()


class SSLEmailBackend(smtp.EmailBackend):
    """
    Required for Amazon SES
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('timeout', 5)
        super(SSLEmailBackend, self).__init__(*args, **kwargs)

    def open(self):
        if self.connection:
            return False
        try:
            logger.debug("sending email via %s" % self.host)
            self.connection = smtplib.SMTP_SSL(self.host, self.port,
                                               local_hostname=DNS_NAME.get_fqdn())
            if self.username and self.password:
                self.connection.login(self.username, self.password)
            return True
        except Exception, exc:
            logger.error(exc)
            if not self.fail_silently:
                raise