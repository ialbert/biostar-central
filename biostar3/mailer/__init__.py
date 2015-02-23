"""
A simplifed version of https://github.com/BradWhittington/django-templated-email

Allows an email to be specified in a template.

Extracts all information from a single template like so. All three blocks must be present.

    {% block subject %} declares the subject
    {% block text %} declares text/plain
    {% block html %} declares text/html

"""

from django.template.loader import get_template
from django.template.loader_tags import BlockNode
from django.core.mail import EmailMultiAlternatives
from django.template import Context
import logging

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
        self.subj = get_node(template, name=self.SUBJECT).render(context)
        self.text = get_node(template, name=self.TEXT).render(context)
        self.html = get_node(template, name=self.HTML).render(context)

    def send(self, from_email, to, cc=None, bcc=None, headers=None):
        msg = EmailMultiAlternatives(
            subject=self.subj, body=self.text, from_email=from_email, to=to,
            cc=cc, bcc=bcc, headers=headers,
        )
        if self.html:
            msg.attach_alternative(self.html, 'text/html')

        msg.send()