"""
User related tests.

These will execute when you run "manage.py test".
"""
from __future__ import print_function, unicode_literals, absolute_import, division
import logging
from django.conf import settings
from biostar.apps.users.models import User, Profile
from django.test import TestCase

logging.disable(logging.INFO)

class UserTest(TestCase):
    def test_user_creation(self):
        """
        Testing users and their profile creation
        """
        eq = self.assertEqual

        # Create a new usr
        user = User.objects.create(email="foo@bar.com")

        # A user will automatically get a profile
        eq (user.profile.user_id, user.id)
