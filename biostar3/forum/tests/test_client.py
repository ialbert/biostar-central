from __future__ import absolute_import, division, unicode_literals

from django.test import TestCase
from django.conf import settings

from biostar3.forum import apps, models
from biostar3.forum.models import User
from django.core.urlresolvers import reverse

from faker import Factory

from django.test import Client


class ClientTests(TestCase):

    def test_navigation(self):
        """
        Tests client navigation.
        """
        host = "www.localhost.com"

        c = Client(HTTP_HOST=host)
        r = c.get('/')

        self.assertEqual(r.status_code, 200)

        r = c.post(reverse("search"), {'q': 'blast'})

        self.assertEqual(r.status_code, 302)


