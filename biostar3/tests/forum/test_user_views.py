from django.core.urlresolvers import reverse

from biostar3.tests.client_test_case import ClientTestCase
from biostar3.forum.models import User

import logging
logging.disable(logging.INFO)

from faker import Factory
f = Factory.create()

class LoginTests(ClientTestCase):
    def test_blank(self):
        response = self.post("sign_up", data=dict(), follow=True)

        self.assertRedirects(response, reverse('account_login'), host=self.HOST)

    def test_no_email(self):
        response = self.post("sign_up", data=dict(password=f.word()))

        self.assertRedirects(response, reverse('account_login'), host=self.HOST)

    def test_no_password(self):
        response = self.post("sign_up", data=dict(email=f.email()))

        self.assertRedirects(response, reverse('account_login'), host=self.HOST)

    def test_valid_user(self):
        user = User.objects.create(name=f.name(), email=f.email())
        password = f.word()
        user.set_password(password)
        user.save()

        response = self.post("sign_up", data=dict(email=user.email, password=password))

        self.assertRedirects(response, reverse('home'), host=self.HOST)

    def test_invalid_user(self):
        response = self.post("sign_up", data=dict(email=f.email(), password=f.word()))

        self.assertRedirects(response, reverse('account_login'), host=self.HOST)

class SignUpTests(ClientTestCase):
    def test_blank(self):
        response = self.post("sign_up", data=dict(signup="1"), follow=True)

        self.assertRedirects(response, reverse('account_login'), host=self.HOST)

    def test_no_email(self):
        response = self.post("sign_up", data=dict(signup="1", password=f.word()))

        self.assertRedirects(response, reverse('account_login'), host=self.HOST)

    def test_no_password(self):
        response = self.post("sign_up", data=dict(signup="1", email=f.email()))

        self.assertRedirects(response, reverse('account_login'), host=self.HOST)

    def test_user_exists(self):
        user = User.objects.create(name=f.name(), email=f.email())
        password = f.word()
        user.set_password(password)
        user.save()

        response = self.post("sign_up", data=dict(signup="1", email=user.email, password=password))

        # Designed behavior is to act as login
        self.assertRedirects(response, reverse('home'), host=self.HOST)

    def test_user_does_not_exist(self):
        email = f.email()
        password = f.word()

        # User does not exist yet
        self.assertEqual(User.objects.filter(email=email).count(), 0)

        response = self.post("sign_up", data=dict(signup="1", email=email, password=password))

        # Successfully login and user exists
        self.assertRedirects(response, reverse('home'), host=self.HOST)
        self.assertEqual(User.objects.filter(email=email).count(), 1)