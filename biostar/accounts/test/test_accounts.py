import logging, os
from django.test import TestCase
from django.test import Client
from biostar.engine import auth
from biostar.engine import models

from django.urls import reverse


logger = logging.getLogger('engine')

class UserAccountTests(TestCase):

    def setUp(self):
        self.user = models.User.objects.filter(is_superuser=True).first()


    def visit_urls(self, urls, code):
        c = Client()
        for url in urls:
            resp = c.get(url)
            if resp.status_code != code:
                # print (resp.content)
                # We already know it is an error.
                # Use this to prints the url and the code.
                logger.error(f"")
                logger.error(f"Error accessing: {url}, code={resp.status_code}")
                self.assertEqual(url, code)

    def test_page_responses(self):

        urls = [
            reverse('index'), reverse('info'), reverse('logout'),
            reverse('login'), reverse('signup'),
            reverse('profile', kwargs=dict(id=self.user.id))
        ]

        self.visit_urls(urls, 200)


class PasswordResetTest(TestCase):


    def setUp(self):
        return


    #def test_password_reset(self):
    #    return







