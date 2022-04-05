import logging, os
from unittest.mock import patch, MagicMock
from django.test import TestCase, override_settings
from django.test import Client
from biostar.accounts import models, views, auth
from django.core import signing

from django.conf import settings
from django.urls import reverse

from biostar.utils.helpers import fake_request, get_uuid


logger = logging.getLogger('engine')


@override_settings(RECAPTCHA_PRIVATE_KEY="", RECAPTCHA_PUBLIC_KEY="")
class UserAccountTests(TestCase):

    def setUp(self):
        self.password = "testing"
        self.user = models.User.objects.create_user(username=f"user1", email="tested@l.com")
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

    def test_user_logout(self):
        "Test user logout interface."

        c = Client()
        c.login(username=self.user.username, email=self.user.email,
                password=self.password)

        url = reverse("logout")
        resp = c.post(url, data={})

        self.assertEqual(resp.status_code, 302)

    def test_page_responses(self):

        urls = [
            reverse("inbox"),
            reverse('accounts_index'),
            reverse('logout'),
            reverse('login'),
            reverse("user_moderate", kwargs=dict(uid=self.user.id)),
            reverse("user_profile", kwargs=dict(uid=self.user.profile.uid)),
            reverse('edit_profile'),
            reverse('password_reset'),
            reverse('password_reset_done'),
            reverse('password_reset_complete'),
        ]

        self.visit_urls(urls, 200)

    def test_redirect(self):
        urls = [
        ]
        self.visit_urls(urls, 302)


@override_settings(RECAPTCHA_PRIVATE_KEY="", RECAPTCHA_PUBLIC_KEY="")
class LoginTest(TestCase):


    def setUp(self):
        self.password = "testing"
        self.user = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com")
        self.user.set_password(self.password)
        self.user.save()

    def test_login(self):
        "Test Valid login"
        data = {"email": self.user.email, "password":self.password}
        url = reverse("login")

        c = Client()
        resp = c.post(url, data=data)

        self.assertEqual(resp.status_code, 302)
        self.assertTrue(resp.url == settings.LOGIN_REDIRECT_URL,
                         f"Invlaid redirection when logging in.\nexpected: /\ngot:{resp.url}")

        return

    def test_invalid_creds(self):
        "Test invalid email and password"
        data1 = {"email": "foo", "password": self.password}
        data2 = {"email": self.user.email, "password": "bar"}
        url = reverse("login")
        for tests in (data1, data2):

            c = Client()
            resp = c.post(url, data=tests)

            self.assertEqual(resp.status_code, 200)


@override_settings(RECAPTCHA_PRIVATE_KEY="", RECAPTCHA_PUBLIC_KEY="")
class SignUpTest(TestCase):


    def setUp(self):
        self.password = "testing"
        self.email = "tested@email.com"

    @override_settings(ALLOW_SIGNUP=True)
    def test_signup(self):
        "Test the signup interface."

        valid = {"email": self.email, "password1":self.password, "password2": self.password, 'code':302}
        invalid = {"email": self.email, "password1": self.password, "password2": "Fail", 'code':200}
        signup_url = reverse("signup")

        for test in (valid, invalid):
            code = test['code']
            del test['code']

            c = Client()
            resp = c.post(signup_url, data=test)
            self.assertEqual(resp.status_code, code)

        # Test logging in after signup
        login_url = reverse("login")
        login_data = {"email": self.email, "password": self.password}
        resp = c.post(login_url, data=login_data)
        self.assertEqual(resp.status_code, 302)

        #views.user_login()


@override_settings(RECAPTCHA_PRIVATE_KEY="", RECAPTCHA_PUBLIC_KEY="")
class ProfileTest(TestCase):


    def setUp(self):
        self.password = "testing"
        self.user = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com")

        self.user.set_password(self.password)
        self.user.save()

        return


    def test_profile(self):
        "Test profile with a logged in user with GET Request"
        data = {}
        url = reverse("user_profile", kwargs=dict(uid=self.user.profile.uid))

        request = fake_request(url=url, data=data, user=self.user, method="GET")

        response = views.user_profile(request=request, uid=self.user.profile.uid)
        self.assertEqual(response.status_code, 200, "Can not load user profile")


    @patch('biostar.accounts.models.User', MagicMock(name="save"))
    def test_edit_profile(self):
        "Test editing profile with POST request."

        data = {"email":"new@new.com", "name":"new name", "handle":"new",
                "digest_prefs" : models.Profile.DAILY_DIGEST,
                "message_prefs" : models.Profile.LOCAL_MESSAGE}

        url = reverse("edit_profile")

        request = fake_request(url=url, data=data, user=self.user)

        response = views.edit_profile(request=request)

        self.assertEqual(response.status_code, 302, "Can not redirect after editing profile")

    def test_notify(self):
        "Test the notification toggling"

        url = reverse("toggle_notify")

        request = fake_request(url=url, data={}, user=self.user)

        response = views.toggle_notify(request=request)

        self.assertEqual(response.status_code, 302)



    def test_banned_user_login(self):
        "Test banned user can not login "

        self.user.profile.state = models.Profile.BANNED
        self.user.profile.save()

        message, valid = auth.validate_login(email=self.user.email, password=self.password)

        print(message, valid)

        self.assertTrue(not valid)
