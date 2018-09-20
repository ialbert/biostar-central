import logging, os
from django.test import TestCase
from django.test import Client
from biostar.forum import auth, models
from biostar.accounts.models import User

from django.urls import reverse


logger = logging.getLogger('engine')


class ForumNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        self.owner = User.objects.create(username="test", email="test@test.com")
        self.owner.set_password("testing")

        # Create a test post
        self.post = auth.create_post(title="Test", author=self.owner, content="Test",
                                     post_type=models.Post.QUESTION)

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
        "Checking public forum pages"

        urls = [
            reverse('post_list'),
            reverse('post_create'),
            reverse("community_list"),
            reverse("tags_list"),

            reverse('post_view', kwargs=dict(uid=self.post.uid)),
            reverse('post_edit', kwargs=dict(uid=self.post.uid)),

        ]

        self.visit_urls(urls, [200])

    def test_page_redirect(self):
        "Testing that a redirect occurs for some pages"
        urls = [
            reverse("post_comment")
        ]

        self.visit_urls(urls, [302])



