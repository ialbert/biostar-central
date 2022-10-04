import logging
import re
import textwrap

from django.core.mail import EmailMultiAlternatives
from django.core.mail import send_mail, send_mass_mail, get_connection
from django.template import Context, Template
from django.template.loader import get_template
from django.conf import settings

logger = logging.getLogger("engine")

# Pattern to extract named blocks from a django template.
block_patt = r'{%\s+block\s+(?P<name>(.+?))\s+%}(?P<value>(.+?)){%\s+endblock\s+%}'
block_regx = re.compile(block_patt)


def strip(text):
    return text.strip()


def get_block(text, name):
    for match in re.finditer(block_patt, text, re.DOTALL | re.MULTILINE | re.IGNORECASE):
        found = match.group("name")
        value = match.group("value")
        if found == name:
            return Template(value)

    msg = "*** Error: missing block name: {}".format(name)
    logger.error(msg)
    return Template(msg)


def safe_render(templ, context):
    try:
        return templ.render(Context(context))
    except Exception as exc:
        msg = f"*** template error: {exc}"
        logger.error(msg)
        return msg


def first_line(text):
    """
    Returns the first non-empty line of a multi-line text
    """
    lines = map(strip, text.splitlines())
    lines = list(filter(None, lines))
    first = lines[0] if lines else ''
    return first


class EmailTemplate(object):
    """
    Generates a subject, text and html based email from a single template.
    """

    def __init__(self, name):
        self.template = get_template(name)
        self.content = open(self.template.origin.name).read()
        self.subj = get_block(self.content, "subject")
        self.text = get_block(self.content, "text")
        self.html = get_block(self.content, 'html')

    def render(self, context):
        subj = safe_render(self.subj, context)
        text = safe_render(self.text, context)
        html = safe_render(self.html, context) if self.html else ''
        subj = first_line(subj)
        return subj, text, html

    def send(self, context, from_email, recipient_list):

        recipients = ", ".join(recipient_list)

        # Skip sending emails during data migration
        if settings.DATA_MIGRATION:
            logger.info(f"skip email to: {recipients} DATA_MIGRATION={settings.DATA_MIGRATION} ")
            return
        else:
            logger.info(f"sending email to: {recipients}")

        subject, text, html = self.render(context)
        # Text may be indented in template.
        text = textwrap.dedent(text)

        if len(html) > 10:
            send_html_mail(
                subject=subject,
                message=text,
                message_html=html,
                from_email=from_email,
                recipient_list=recipient_list)
        else:
            send_mail(
                subject=subject,
                message=text,
                from_email=from_email,
                recipient_list=recipient_list,
                html_message=html)

    def send_mass(self, context, from_email, recipient_list):
        """
        Send mass individual mail to list of recipients
        """
        subject, text, html = self.render(context)

        # Text may be indented in template.
        text = textwrap.dedent(text)

        # Send mass html email
        if len(html) < 10:
            # Format mass mail
            datatuple = ((subject, text, from_email, [rec]) for rec in recipient_list)
            send_mass_mail(datatuple=datatuple, fail_silently=False)
        else:
            send_mass_html_mail(subject=subject,
                                message=text,
                                message_html=html,
                                from_email=from_email,
                                recipient_list=recipient_list)


def send_mass_html_mail(subject, message, message_html, from_email, recipient_list):
    """

    """
    connection = get_connection(fail_silently=False)

    def make_email(rec):
        msg = EmailMultiAlternatives(subject=subject,
                                     body=message,
                                     from_email=from_email,
                                     to=[rec],
                                     connection=connection)
        msg.attach_alternative(message_html, "text/html")
        return msg

    # Format mass mail
    messages = list(map(make_email, recipient_list))

    return connection.send_messages(messages)


def send_html_mail(subject, message, message_html, from_email, recipient_list):
    """
    Sends an HTML email.
    """
    msg = EmailMultiAlternatives(subject, message, from_email, recipient_list)
    msg.attach_alternative(message_html, "text/html")
    msg.send(fail_silently=False)
