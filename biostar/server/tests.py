from django.test import TestCase, SimpleTestCase
from django.test import Client
from django.core.urlresolvers import reverse
from django.conf import settings
from biostar.apps.users.models import User

import logging

logging.disable(logging.WARNING)

class UserTest(SimpleTestCase):

    # The name of test users
    NAME_1, EMAIL_1, PASSWD_1 = "John Doe", "john@example.org", "0123567"
    NAME_2, EMAIL_2, PASSWD_2 = "Jane Doe", "jane@example.org", "3456789"

    # The name of test posts
    TITLE_1, CAT_1, TAG_VAL_1 = "Post 1", "Job", "tagA tagB galaXY"
    TITLE_2, CAT_2, TAG_VAL_2 = "Post 2", "galAxy", "tagA tagB tagD"

    def code(self, response, code=200):
        self.assertEqual(response.status_code, code)

    def setUp(self):
        # Sign up a user
        r = self.client.get(reverse("signup"), follow=True)

        # Syncdb add the admin user.
        self.assertEqual(User.objects.all().count(), 1)

        # Sign up user 1.
        r = self.client.post(reverse("account_signup"), dict(email=self.EMAIL_1, password1=self.PASSWD_1, password2=self.PASSWD_1), follow=True)
        self.assertContains(r, "My Tags")
        self.assertEqual(User.objects.all().count(), 2)

        # Logout user 1.
        r = self.client.get(reverse("logout"), follow=True)
        self.assertNotContains(r, "My Tags")

        # Sign up user 2.
        r = self.client.post(reverse("account_signup"), dict(email=self.EMAIL_2, password1=self.PASSWD_2, password2=self.PASSWD_2), follow=True)
        self.assertContains(r, "My Tags")
        self.assertEqual(User.objects.all().count(), 3)

         # Logout user 2.
        r = self.client.get(reverse("logout"), follow=True)
        self.assertNotContains(r, "My Tags")

    def test_user_login(self):
        eq = self.assertEqual


        r = self.client.post(reverse("account_login"), dict(login=self.EMAIL_1, password=self.PASSWD_1), follow=True)
        self.assertContains(r, "My Tags")


        # Go to the
        #r = self.client.get(reverse("home"))
        #self.assertContains(r, "My Tags")


class SiteTest(TestCase):

    def code(self, response, code=200):
        self.assertEqual(response.status_code, code)

    def test_site_navigation(self):
        "Testing site navigation"
        eq = self.assertEqual

        # Main site navigation.
        names = "home user-list tag-list help about faq policy rss latest-feed".split()
        for name in names:
            r = self.client.get(reverse(name))
            self.code(r)

        # Check that default categories work.
        for topic in settings.CATEGORIES:
            r = self.client.get(reverse("topic-list", kwargs=dict(topic=topic)))
            self.code(r)

    def test_redirects(self):
        "Testing page redirects"

        # Pages with redirects.
        names = "login logout signup new-post".split()
        for name in names:
            r = self.client.get(reverse(name))
            self.code(r, 302)

        # Pages that take parameters and redirect.
        names = "user-edit post-edit".split()
        for name in names:
            r = self.client.get(reverse(name, kwargs=dict(pk=1)))
            self.code(r, 302)



