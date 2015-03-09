from __future__ import absolute_import, division, unicode_literals

from django.test import TestCase
from django.conf import settings

from biostar3.forum import apps, models
from biostar3.forum.models import User
from django.core.urlresolvers import reverse

from faker import Factory

from django.test import Client

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
            self.assertTrue(re.search(pattern, r.content, re.IGNORECASE))

        return r


    def test_navigation(self):
        """
        Tests client navigation.
        """

        c = Client(HTTP_HOST=HOST)
        r = self.get(c, "home")

        r = self.post(c, "search", data={'q': 'blast'})


    def test_user_signup(self):
        EQ = self.assertEqual

        c = Client(HTTP_HOST=HOST)
        r = self.get(c, "account_login", pattern='social authentication')

        data = dict(email=ADMIN_EMAIL, password=settings.SECRET_KEY)

        r = self.post(c, "sign_up", data=data, follow=True, pattern="success")






