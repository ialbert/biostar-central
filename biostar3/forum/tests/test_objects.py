from __future__ import absolute_import, division, print_function, unicode_literals

from django.test import TestCase
from django.conf import settings

from biostar3.forum import apps
from biostar3.forum.models import User

class SimpleTests(TestCase):

    def test_admin_user(self):
        """
        Admin user must exists and has the right permissions.
        """
        for user, email in settings.ADMINS:
            user = User.objects.get(email=email)
            funcs = [
                apps.authorize_post_mod, apps.authorize_user_ban, apps.authorize_user_mod
            ]
            for func in funcs:
                self.assertEqual(func(user), True)