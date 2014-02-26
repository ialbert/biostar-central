from django.test import TestCase
from django.test import Client
from django.core.urlresolvers import reverse
from django.conf import settings
from biostar.apps.users.models import User

import logging

logging.disable(logging.WARNING)

class UserTest(TestCase):

    def setUp(self):
        # Sign up a user
        c = Client()
        r = c.get(reverse("signup"), follow=True)


    def test_user_login(self):
        pass


class SiteTest(TestCase):

    def code(self, response, code=200):
        self.assertEqual(response.status_code, code)

    def test_site_navigation(self):
        """
        Testing site navigation
        """
        eq = self.assertEqual
        c = Client()

        # Main site navigation.
        names = "home user-list tag-list help about faq policy rss latest-feed".split()
        for name in names:
            r = c.get(reverse(name))
            self.code(r)

        # Pages with redirects.
        names = "login logout signup new-post".split()
        for name in names:
            r = c.get(reverse(name))
            self.code(r, 302)

        # Pages that take parameters and redirect.
        names = "user-edit post-edit".split()
        for name in names:
            r = c.get(reverse(name, kwargs=dict(pk=1)))
            self.code(r, 302)

        # Check that default categories work.
        for topic in settings.CATEGORIES:
            r = c.get(reverse("topic-list", kwargs=dict(topic=topic)))
            self.code(r)

