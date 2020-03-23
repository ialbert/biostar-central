import logging

from django.test import TestCase, override_settings
from django.shortcuts import reverse
from biostar.accounts.models import User
from biostar.utils.helpers import fake_request
from biostar.accounts import views, forms

@override_settings(RECAPTCHA_PRIVATE_KEY="", RECAPTCHA_PUBLIC_KEY="")
class ModerateUser(TestCase):
    def setUp(self):

        self.user1 = User.objects.create(username=f"foo1", email="foo@tested.com",
                                         password="foo", is_superuser=True, is_staff=True)
        self.user2 = User.objects.create(username=f"foo2", email="foo2@tested.com",
                                       password="foo2")
        pass

    def test_user_moderate(self):
        "Test user moderation"

        for action, _ in forms.UserModerate.CHOICES:

            url = reverse("user_moderate", kwargs=dict(uid=self.user2.id))

            data = {"action": action}
            request = fake_request(user=self.user1, data=data, url=url)

            response = views.user_moderate(request=request, uid=self.user2.id)

            self.assertTrue(response.status_code == 302, "Error moderating user.")

        pass

    def test_debug_user(self):
        "Test logging in as a user"

        url = reverse("debug_user")

        request = fake_request(user=self.user1, data=dict(uid=self.user2.profile.uid), url=url)

        response = views.debug_user(request=request)

        self.assertTrue(response.status_code == 302, "Error debugging user.")
