"""
This file demonstrates writing tests using the unittest module. These will pass
when you run "manage.py test".

Replace this with more appropriate tests for your application.
"""
import logging
from django.conf import settings
from biostar.apps.accounts.models import User, Profile
from django.test import TestCase

logging.disable(logging.CRITICAL)

class SimpleTest(TestCase):
    def test_users(self):
        """
        Testing users and their profile creation
        """
        eq = self.assertEqual
        admin = User.objects.get(email=settings.ADMIN_EMAIL)

        # An admin user is created by default.
        eq(admin.name, settings.ADMIN_NAME)

        # It must have a profile created.
        eq(admin.profile.user.email, settings.ADMIN_EMAIL)
