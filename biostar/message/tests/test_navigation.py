import logging, os
from django.test import TestCase
from django.test import Client
from biostar.message import auth, models
from biostar.accounts.models import User, Profile

from django.urls import reverse
import uuid


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


logger = logging.getLogger('engine')


class ForumNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        self.owner = User.objects.create(username=f"tested{get_uuid(10)}", email="tested@tested.com")
        self.owner.set_password("tested")

        # Create a tested message
        body = "Test message"
        rec_list = [self.owner]
        messages = auth.create_local_messages(body=body, sender=self.owner, rec_list=rec_list, subject="",
                                              parent=None, source=models.Message.REGULAR,
                                              mtype=Profile.LOCAL_MESSAGE)
        self.message = messages[0]
        self.owner.save()

    def visit_urls(self, urls, codes):
        c = Client()
        c.login(username=self.owner.username, email='tested@tested.com', password='tested')
        for url in urls:
            print(url)
            resp = c.get(url)
            code = resp.status_code
            if code not in codes:
                # We already know it is an error.
                # Use this to prints the url and the code.
                logger.error(f"Error accessing: {url}, code={code} not in expected values {codes}")
                self.assertTrue(code in codes)

    def test_public_pages(self):
        "Checking messages pages"

        urls = [
            reverse('inbox'),
            reverse('outbox'),
            reverse('compose'),
            reverse("inbox_view", kwargs=dict(uid=self.message.uid)),
            reverse('outbox_view', kwargs=dict(uid=self.message.uid)),


        ]

        self.visit_urls(urls, [200])