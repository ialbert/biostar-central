import logging

from django.test import TestCase
from django.urls import reverse

from biostar.accounts.models import User

from biostar.forum import models, views, auth, forms, const
from biostar.forum.tests.util import fake_request
from biostar.forum.util import get_uuid

logger = logging.getLogger('engine')


class PostTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"tested{get_uuid(10)}", email="tested@tested.com",
                                         password="tested", is_superuser=True, is_staff=True)
        self.owner.save()

    def test_ajax_create(self):
        return

    def test_ajax_edit(self):
        return

    def test_ajax_search(self):
        return

    def test_inplace_form(self):
        return

    def test_ajax(self):
        return

