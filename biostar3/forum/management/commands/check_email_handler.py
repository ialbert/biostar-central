from __future__ import print_function, unicode_literals, absolute_import, division
import logging, requests, pyzmail

from django.core.management.base import BaseCommand
from django.conf import settings
from django.core.urlresolvers import reverse

from biostar3.forum import auth, models


logger = logging.getLogger('biostar')

from faker import Factory

faker = Factory.create()

def compose_new(user):

    target = settings.EMAIL_ADDRESS_PATTERN.format("new", user.profile.uuid, settings.DEFAULT_GROUP_DOMAIN)

    sender = (u'Me', 'me@foo.com')
    recipients = [ target ]
    subject = faker.bs()
    text_content = u'Bonjour aux *Fran\xe7ais*!'
    prefered_encoding = 'iso-8859-1'
    text_encoding = 'iso-8859-1'
    html = None

    payload, mail_from, rcpt_to, msg_id = pyzmail.compose_mail(
        sender, recipients, subject, prefered_encoding,
        (text_content, text_encoding), html=html
    )

    return payload

class Command(BaseCommand):
    help = 'Checks the email handler'

    def handle(self, *args, **options):
        # Get the admin user
        site = models.Site.objects.get_current()

        admin = models.User.objects.filter(email=settings.ADMINS[0][1]).first()

        email_handler = "http://{}{}".format(site.domain, reverse("email_handler"))

        models.Post.objects.filter(root=None).delete()

        content = compose_new(admin)

        print(content)

        data = dict(key=settings.EMAIL_HANDLER_SECRET_KEY, content=content)

        r = requests.post(email_handler, data=data)

        print('-' * 10)
        print(r.text)

