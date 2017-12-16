import logging, os
from unittest.mock import patch, MagicMock
from django.test import TestCase
from django.test import Client
from biostar.accounts import models, views

from django.urls import reverse

from . import util


logger = logging.getLogger('engine')

class UserAccountTests(TestCase):

    def setUp(self):
        self.password = "testing"
        self.user = models.User.objects.create_user(username="test", email="test@l.com")
        self.user.set_password(self.password)
        self.user.save()

    def visit_urls(self, urls, code):

        c = Client()
        c.login(username=self.user.username, email=self.user.email,
                password=self.password)

        for url in urls:
            resp = c.get(url)
            if resp.status_code != code:
                # print (resp.content)
                # We already know it is an error.
                # Use this to prints the url and the code.
                logger.error(f"Error accessing: {url}, code={resp.status_code}")
                self.assertEqual(url, code)

    def test_page_responses(self):

        urls = [
            reverse('index'),
            reverse('logout'),
            reverse('login'),
            reverse('profile'),
            reverse('edit_profile'),
            reverse('password_reset'),
            reverse('password_reset_done'),
            reverse('password_reset_complete'),
        ]

        self.visit_urls(urls, 200)

    def test_redirect(self):
        urls = [

            reverse('signup')
        ]
        self.visit_urls(urls, 302)


class LoginTest(TestCase):


    def setUp(self):
        self.password = "testing"
        self.user = models.User.objects.create_user(username="test", email="test@l.com")
        self.user.set_password(self.password)

    def test_login(self):
        data = {"email": self.user.email, "password":self.password}
        url = reverse("login")

        request = util.fake_request(url=url, data=data, user=self.user)

        response = views.user_login(request=request)

        self.assertEqual(response.status_code, 302)

        return


class ProfileTest(TestCase):


    def setUp(self):
        self.password = "testing"
        self.user = models.User.objects.create_user(username="test", email="test@l.com")
        self.user.set_password(self.password)

        return


    def test_profile(self):
        "Test profile with a logged in user with GET Request"
        data = {}
        url = reverse("profile")

        request = util.fake_request(url=url, data=data, user=self.user, method="GET")

        response = views.profile(request=request)
        self.assertEqual(response.status_code, 200, "Can not load user profile")


    @patch('biostar.accounts.models.User', MagicMock(name="save"))
    def test_edit_profile(self):
        "Test editing profile with POST request"

        data = {"email":"new@new.com", "first_name":"new name"}

        url = reverse("edit_profile")

        request = util.fake_request(url=url, data=data, user=self.user)

        response = views.edit_profile(request=request)

        self.assertEqual(response.status_code, 302, "Can not redirect after editing profile")
        self.assertTrue(reverse("profile") == response.url,
                        "Could not redirect to correct page after editing")
