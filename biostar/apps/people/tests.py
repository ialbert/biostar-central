"""
User related tests.

These will execute when you run "manage.py test".
"""
import logging
from django.conf import settings
from biostar.apps.people.models import User, Profile
from django.test import TestCase

logging.disable(logging.INFO)

class UserTest(TestCase):
    def test_user_creation(self):
        """
        Testing users and their profile creation
        """
        eq = self.assertEqual
        admin = User.objects.get(email=settings.ADMIN_EMAIL)

        # An admin user is created by default.
        eq(admin.name, settings.ADMIN_NAME)

        # It must have a profile created.
        eq(admin.profile.user.email, settings.ADMIN_EMAIL)
