import logging, os
from unittest.mock import patch, MagicMock
from django.test import TestCase
from django.test import Client
from biostar.accounts import models, views, auth

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
        self.user.save()

    def test_login(self):
        "Test Valid login"
        data = {"email": self.user.email, "password":self.password}
        url = reverse("login")

        c = Client()
        resp = c.post(url, data=data)

        self.assertEqual(resp.status_code, 302)
        self.assertTrue(resp.url == "/",
                         f"Invlaid redirection when logging in.\nexpected: /\ngot:{resp.url}")

        return

    def test_invalid_creds(self):
        "Test invalid email and password"
        data1 = {"email": "foo", "password": self.password}
        data2 = {"email": self.user.email, "password": "bar"}

        for tests in (data1, data2):

            url = reverse("login")

            c = Client()
            resp = c.post(url, data=tests)

            self.assertEqual(resp.status_code, 200)



class ProfileTest(TestCase):


    def setUp(self):
        self.password = "testing"
        self.user = models.User.objects.create_user(username="test", email="test@l.com")

        self.user.set_password(self.password)
        self.user.save()

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

        data = {"email":"new@new.com", "name":"new name"}

        url = reverse("edit_profile")

        request = util.fake_request(url=url, data=data, user=self.user)

        response = views.edit_profile(request=request)

        self.assertEqual(response.status_code, 302, "Can not redirect after editing profile")


    def test_notify(self):
        "Test the notification toggling"

        url = reverse("toggle_notify")

        request = util.fake_request(url=url, data={}, user=self.user)

        response = views.toggle_notify(request=request)

        self.assertEqual(response.status_code, 302)


    def test_banned_user_login(self):
        "Test banned user can not login "

        self.user.profile.state = models.Profile.BANNED
        self.user.profile.save()

        message, valid = auth.check_user(email=self.user.email, password=self.password )

        print(message, valid)

        self.assertTrue(not valid)








