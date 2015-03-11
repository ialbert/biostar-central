from __future__ import absolute_import, division, unicode_literals

from django.test import TestCase
from django.conf import settings

from biostar3.forum import apps, models
from biostar3.forum.models import Post

from biostar3.forum.models import User
from django.core.urlresolvers import reverse

from faker import Factory

from django.test import Client
import random, faker
HOST = "www.localhost.com"

ADMIN_NAME, ADMIN_EMAIL = settings.ADMINS[0]

import logging, re

logging.disable(logging.ERROR)


class ClientTests(TestCase):
    def get(self, client, url, code=200, pattern=None):
        r = client.get(reverse(url))
        self.assertEqual(r.status_code, code)
        if pattern:
            self.assertTrue(re.search(pattern, r.content, re.IGNORECASE))
        return r

    def post(self, client, url, data, pattern=None, follow=False):
        r = client.post(reverse(url), data, follow=follow)

        if follow:
            self.assertEqual(r.status_code, 200)
        else:
            self.assertEqual(r.status_code, 302)

        if pattern:
            result = re.search(pattern, r.content, re.IGNORECASE)
            if not result:
                print "Unable to find %s pattern in content." % (pattern)
                self.assertTrue(result)

        return r

    def test_navigation(self):
        """
        Tests client navigation.
        """

        c = Client(HTTP_HOST=HOST)
        r = self.get(c, "home")
        r = self.post(c, "search", data={'q': 'blast'}, follow=True)


    def test_user_posting(self):
        """
        Test user signup
        """
        EQ = self.assertEqual

        c = Client(HTTP_HOST=HOST)
        r = self.get(c, "account_login", pattern='social authentication')

        f = faker.Factory.create()
        email =  f.email()

        data = dict(email=email , password=email, signup="1")

        r = self.post(c, "sign_up", data=data, follow=True, pattern="success")

        r = self.get(c, "new_post")

        title = f.sentence()[:100]
        content = f.text()
        tags = "hello, world"

        post_type = models.Post.QUESTION

        data = dict(title=title , tags=tags, content=content, type=post_type)

        r = self.post(c, "new_post", data=data, follow=True, pattern=title)

        self.assertTrue(Post.objects.filter(title=title))



