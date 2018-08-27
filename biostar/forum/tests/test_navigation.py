import logging, os
from django.test import TestCase
from django.test import Client
from biostar.engine import auth
from biostar.engine import models

from django.urls import reverse


logger = logging.getLogger('engine')


class ForumNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        self.owner = models.User.objects.create(username="test", email="test@test.com")
        self.owner.set_password("testing")
        self.owner.save()

    def visit_urls(self, urls, codes):
        c = Client()
        c.login(username="test", email='test@test.com', password='testing')
        for url in urls:
            resp = c.get(url, data={"q":"test"})
            code = resp.status_code
            if code not in codes:
                # We already know it is an error.
                # Use this to prints the url and the code.
                logger.error(f"")
                logger.error(f"Error accessing: {url}, code={code} not in expected values {codes}")
                self.assertTrue(code in codes)

    def test_public_pages(self):
        "Checking public pages"

        urls = [
            reverse('post_list'),
            # reverse('post_view', kwargs=dict()),
            # reverse('post_create', kwargs=dict()),
            # reverse("subs_action", kwargs=dict()),
            # reverse("post_comment", kwargs=dict()),
            # reverse("update_vote", kwargs=dict()),
            # reverse("tags_list", kwargs=dict()),
            # reverse("community_list", kwargs=dict())


        ]

        self.visit_urls(urls, [200])

    def test_page_redirect(self):
        "Testing that a redirect occurs for some pages"
        urls = [

        ]

        self.visit_urls(urls, [302, 200])



