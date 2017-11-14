import logging, os
from django.test import Client
from django.test import TestCase
from biostar.engine import models
from biostar.engine import auth

logger = logging.getLogger('engine')


class AccountsNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        user = models.User.objects.all().first()
        #self.project = auth.create_project(user=user, name="Test project")

