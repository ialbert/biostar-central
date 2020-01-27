import logging, os
from django.test import TestCase
from django.test import Client
from django.core import management
from biostar.forum import auth, models
from biostar.accounts.models import User
from biostar.forum.models import Badge

from django.urls import reverse
import uuid


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


logger = logging.getLogger('engine')


class PlanetNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        self.owner = User.objects.create(username=f"tested{get_uuid(10)}", email="tested@tested.com")
        self.owner.set_password("tested")
        self.badge = Badge.objects.first()
        # Create a tested post
        self.post = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                     type=models.Post.QUESTION)
        management.call_command('populate')

        self.owner.save()

    def visit_urls(self, urls, codes):
        c = Client()
        c.login(username=self.owner.username, email='tested@tested.com', password='tested')
        for url in urls:
            print(url)
            resp = c.get(url, data={"q":"tested"})
            code = resp.status_code
            if code not in codes:
                # We already know it is an error.
                # Use this to prints the url and the code.
                logger.error(f"")
                logger.error(f"Error accessing: {url}, code={code} not in expected values {codes}")
                self.assertTrue(code in codes)

    def test_public_pages(self):
        "Checking public forum pages"

        urls = [
            reverse("blog_list"),

        ]

        self.visit_urls(urls, [200])



