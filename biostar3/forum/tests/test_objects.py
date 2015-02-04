from __future__ import absolute_import, division, unicode_literals

from django.test import TestCase
from django.conf import settings

from biostar3.forum import apps, models
from biostar3.forum.models import User

from faker import Factory


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
        f = Factory.create()
        count = 10
        for i in range(count):
            user = User.objects.create(name=f.name(), email=f.email())
            for func in AUTH_FUNCS:
                self.assertEqual(func(user), False)
            # Create a few groups.
            apps.create_group(name=f.user_name(), user=user)

        self.assertTrue(models.Group.objects.all().count() > count)