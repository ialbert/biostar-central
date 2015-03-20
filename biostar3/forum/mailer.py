"""
A simplifed version of https://github.com/BradWhittington/django-templated-email

Allows an email to be specified in a template.

Extracts all information from a single template like so. All three blocks must be present.

    {% block subject %} declares the subject
    {% block text %} declares text/plain
    {% block html %} declares text/html

"""
from __future__ import absolute_import, division, print_function, unicode_literals

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
from .models import Message, MessageBody, right_now

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
        node = get_node(template, name=name)
        if node:
            return node.render(context)
        else:
            return ""
    except Exception as exc:
        logger.error(exc)
        result = "internal error, render_node=%s" % exc
        return result

class EmailTemplate(object):
    """

    """
    MESG, SUBJECT, TEXT, HTML = "message", "subject", "text", "html"

    def __init__(self, template_name, data={}):

        # Get the base template.
        self.template = get_template(template_name)

        # Extract the blocks for each part of the email.
        subj = render_node(self.template, name=self.SUBJECT, data=data)
        text = render_node(self.template, name=self.TEXT, data=data)
        html = render_node(self.template, name=self.HTML, data=data)
        mesg = render_node(self.template, name=self.MESG, data=data)

        # Email subject may not contain newlines
        lines = map(string.strip, subj.splitlines())
        self.subj = "".join(lines)

        # Text node may be indented in templates. Remove common whitespace.
        self.text = textwrap.dedent(text)
        self.html = html
        self.mesg = mesg

    def create_messages(self, author, targets=[]):
        """
        Creates all the local messages based on the template.
        """

        if not targets:
            # There is nobody to send the message to.
            return

        if not self.mesg:
            logger.error("message body empty")
            return

        # This is the main message body that will be shown for each message.
        body = MessageBody.objects.create(
            author=author, subject=self.subj, content=self.mesg, html=self.html,
        )

        # This generates the messages for the targets
        def message_generator(body):
            now = right_now()
            for target in targets:
                yield Message(user=target.user, body=body, date=now)

        # Bulk insert for all messages
        Message.objects.bulk_create(message_generator(body), batch_size=100)

    def send_email(self, to, from_email=None, cc=None, bcc=None, headers={}, token=None):
        """
        Sends the emails corresponding to this message.
        """

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
        except Exception as exc:
            logger.error(exc)
            if not self.fail_silently:
                raise