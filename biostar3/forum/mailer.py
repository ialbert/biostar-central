"""
A simplifed version of https://github.com/BradWhittington/django-templated-email

Allows an email to be specified in a template.

Extracts all information from a single template like so. All three blocks must be present.

    {% block subject %} declares the subject
    {% block text %} declares text/plain
    {% block html %} declares text/html

"""
import logging, smtplib, string, textwrap

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
    # Missing nodes are ok.
    return ""

def render_node(template, name, data):
    """
    Renders a node of a template
    """
    try:
        context = Context(data, autoescape=False)
        result = get_node(template, name=name).render(context)
    except Exception, exc:
        logger.error(exc)
        result = "internal error, render_node=%s" % exc
    return result

class EmailTemplate(object):
    """

    """
    SUBJECT, TEXT, HTML = "subject", "text", "html"

    def __init__(self, template_name, data={}):

        # Get the base template.
        self.template = get_template(template_name)

        # Extract the blocks for each part of the email.
        subj = render_node(self.template, name=self.SUBJECT, data=data)
        text = render_node(self.template, name=self.TEXT, data=data)
        html = render_node(self.template, name=self.HTML, data=data)

        # Email subject may not contain newlines
        lines = map(string.strip, subj.splitlines())
        self.subj = "".join(lines)

        # Text node may be indented in templates. Remove common whitespace.
        self.text = textwrap.dedent(text)
        self.html = html

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

        # Create the email.
        msg = EmailMultiAlternatives(
            subject=self.subj, body=self.text, from_email=from_email, to=to,
            cc=cc, bcc=bcc, headers=headers,
        )

        # Attach html if exists.
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