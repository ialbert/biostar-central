from __future__ import absolute_import, division, unicode_literals

import logging

from django.test import TestCase
from django.conf import settings
from django.core import mail

from biostar3.forum import apps, models
from biostar3.forum.models import User, UserGroup, GroupPerm, Post

from biostar3.forum import html

from faker import Factory
from django.contrib.sites.models import Site

logging.disable(logging.INFO)

class SimpleTests(TestCase):
    def test_admin_user(self):
        """
        Test admin user creation.
        """
        TRUE = self.assertTrue

        TRUE(len(settings.ADMINS) > 0)
        for user, email in settings.ADMINS:
            user = User.objects.get(email=email)
            TRUE(user.is_admin)

    def test_user_gen(self):
        """
        Tests user creation.
        """
        EQ = self.assertEqual

        f = Factory.create()
        count = 10

        start = UserGroup.objects.all().count()

        uuid = lambda x: models.make_uuid(x)
        for i in range(count):
            user = User.objects.create(name=f.name(), email=f.email())
            # Create groups
            group = UserGroup.objects.create(domain=uuid(8), name=uuid(8), owner=user)

            # Each user now belongs to the default and group and
            # the group they have created.
            self.assertEqual(len(user.groupsub_set.all()), 2)


        self.assertTrue(UserGroup.objects.all().count() == start + count)


        # Test sending automated emails.
        EQ(len(mail.outbox), count)

    def test_local_links(self):
        "Links to local content are reformatted"
        EQ = self.assertEqual
        f = Factory.create()
        site = Site.objects.get_current()
        user = User.objects.create(email=f.email(), name=f.name())
        link = "http://%s/u/%s" % (site.domain, user.id)
        text = "ABC %s ABC" % link
        result = html.sanitize(text, user).strip()
        expect = '<p>ABC <a href="%s">%s</a> ABC</p>' % (link, user.name)
        EQ(expect, result)

    def test_embed_links(self):
        "Links to gist and youtube are embedded"
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

