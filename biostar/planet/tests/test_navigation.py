import logging, os
from django.test import TestCase, override_settings
from django.test import Client
from django.core import management
from biostar.forum import auth, models
from biostar.accounts.models import User
from biostar.forum.models import Badge

from django.urls import reverse
import uuid


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


__MODULE_DIR = os.path.dirname(auth.__file__)
TEST_ROOT = os.path.join(__MODULE_DIR, 'tests')

logger = logging.getLogger('engine')

PLANET_DIR = os.path.abspath(os.path.join(TEST_ROOT, "feeds"))


@override_settings(PLANET_DIR=PLANET_DIR, INIT_PLANET = False)
class PlanetNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        self.owner = User.objects.create(username=f"tested{get_uuid(10)}", email="tested@tested.com")
        self.owner.set_password("tested")
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



