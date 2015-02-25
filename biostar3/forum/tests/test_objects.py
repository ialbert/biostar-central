from __future__ import absolute_import, division, unicode_literals

from django.test import TestCase
from django.conf import settings
from django.core import mail

from biostar3.forum import apps, models
from biostar3.forum.models import User
from biostar3.forum import html

from faker import Factory
from django.contrib.sites.models import Site

AUTH_FUNCS = [
    apps.authorize_post_mod, apps.authorize_user_ban, apps.authorize_user_mod
]


class SimpleTests(TestCase):
    def test_admin_user(self):
        """
        Test admin user creation.
        """
        for user, email in settings.ADMINS:
            user = User.objects.get(email=email)

            for func in AUTH_FUNCS:
                self.assertEqual(func(user), True)

    def test_user_gen(self):
        """
        Test normal user creation
        """
        EQ = self.assertEqual

        f = Factory.create()
        count = 10

        start_group_count = models.Group.objects.all().count()

        for i in range(count):
            user = User.objects.create(name=f.name(), email=f.email())
            for func in AUTH_FUNCS:
                self.assertEqual(func(user), False)
            # Create a few groups.
            models.get_or_create_group(name=f.user_name(), user=user)

        self.assertTrue(models.Group.objects.all().count() == start_group_count + count)
        # Test sending automated emails.
        EQ(len(mail.outbox), count)

    def test_local_links(self):
        EQ = self.assertEqual
        f = Factory.create()
        site = Site.objects.get_current()
        user = User.objects.create(email=f.email(), name=f.name())
        link = "http://%s/u/%s" % (site.domain, user.id)
        text = "ABC %s ABC" % link
        result = html.sanitize(text, user).strip()
        expect = '<p>ABC <a href="%s">%s</a> ABC</p>' % (link, user.name)
        print result
        EQ(expect, result)

    def test_embed_links(self):
        EQ = self.assertEqual
        pairs = [
            # input, expected
            ("abcd", "abcd"),
            ("ABC https://gist.github.com/123 ABC", "ABC %s ABC" % html.get_embedded_gist(123)),
            ("ABC https://www.youtube.com/watch?v=123 ABC", "ABC %s ABC" % html.get_embedded_youtube(123)),
        ]
        for text, expect in pairs:
            result = html.embed_links(text)
            EQ(expect, result)

