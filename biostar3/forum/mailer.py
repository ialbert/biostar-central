"""
A simplifed version of https://github.com/BradWhittington/django-templated-email

Allows an email to be specified in a template.

Extracts all information from a single template like so. All three blocks must be present.

    {% block subject %} declares the subject
    {% block text %} declares text/plain
    {% block html %} declares text/html

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import smtplib
import textwrap

from django.conf import settings
from django.template.loader import get_template
from django.template.loader_tags import BlockNode
from django.core.mail import EmailMultiAlternatives
from django.template import Context
from django.core.mail.utils import DNS_NAME
from django.core.mail.backends import smtp

from .models import *
from . import auth
from biostar3.utils.compat import *

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

def post_notifications(post, local_targets=[], email_targets=[]):
    """
    Generates notifications on a post.
    """
    # Used in the template context.
    site = Site.objects.get_current()

    # Full url to the post.
    post_url = "%s://%s%s" % (settings.SITE_SCHEME, site.domain,
                              reverse("post_view", kwargs=dict(pk=post.id)))

    # Full url to the author of the post.
    user_url = "%s://%s%s" % (settings.SITE_SCHEME, site.domain,
                              reverse("user_view", kwargs=dict(pk=post.author.id)))

    # The context that will be passed to the template.
    context = dict(post=post, site=site, scheme=settings.SITE_SCHEME,
                   post_url=post_url, user_url=user_url,
                   slug=site.domain)

    # This is the body of the message that gets created.
    em = EmailTemplate("post_created_message.html", data=context)

    # Create local messages to all targets.
    em.create_messages(author=post.author, post=post, users=local_targets)

    # Generate the email messages. Will be bulk inserted.
    def token_generator(obj):
        now = right_now()

        # Generate an email for all candidates.
        for user in email_targets:
            token = auth.make_uuid(size=8)
            em.send_email(to=[user.email], token=token)
            yield ReplyToken(user=user, post=obj, token=token, date=now)

    # Insert the reply tokens into the database.
    ReplyToken.objects.bulk_create(token_generator(post), batch_size=100)


class EmailTemplate(object):
    """
    Represents a Django Template that can be sent as a local message or as an email.
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
        lines = map(strip, subj.splitlines())
        self.subj = "".join(lines)

        # Text node may be indented in templates. Remove common whitespace.
        self.text = textwrap.dedent(text)
        self.html = html
        self.mesg = mesg

    def create_messages(self, author, post=None, users=[]):
        """
        Creates all the local messages based on the template.
        """

        if not users:
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
            for user in users:
                yield Message(user=user, body=body, post=post, date=now)

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
            reply_to = settings.EMAIL_ADDRESS_PATTERN.format("reply", token, "submit")
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