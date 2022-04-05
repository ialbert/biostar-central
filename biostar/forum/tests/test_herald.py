import logging
import os
import shutil
import datetime
from django.core import management
from django.urls import reverse
from django.test import TestCase, override_settings
from django.conf import settings
from biostar.forum import models, api, herald, ajax
from biostar.utils.helpers import fake_request
from biostar.accounts.models import User

logger = logging.getLogger('engine')


class HeraldTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"test", email="tested@tested.com", password="tested")
        self.staff_user = User.objects.create(username=f"test2", is_superuser=True, is_staff=True,
                                              email="tested@staff.com", password="tested")

        # Create one accepted link.
        self.link = models.SharedLink.objects.create(author=self.owner, text='test', url='http://bing.com',
                                                     status=models.SharedLink.ACCEPTED)

    def test_submit(self):
        """Test herald submissions"""

        # Create fake request
        data = {'url': "https://google.com", 'text': 'This is a test link'}

        request = fake_request(url=reverse('herald_list'), data=data, user=self.owner)
        response = herald.herald_list(request=request)

        self.assertEqual(response.status_code, 302, f"Could not redirect :\nresponse:{response}")

    def test_publish(self):
        """ Test herald publish """

        # Create fake request
        request = fake_request(url=reverse('herald_publish'), data={}, method='GET', user=self.staff_user)
        response = herald.herald_publish(request=request)

        self.assertNotEqual(response.url, reverse('herald_list'), "could not publish ")

        pass

    def test_reviews(self):
        """ Test accept/decline workflow for herald."""

        # Create fake request
        data = dict(status='accept')
        request = fake_request(url=reverse('herald_update', kwargs=dict(pk=self.link.pk)), data=data,
                               user=self.staff_user)
        response = ajax.herald_update(request=request, pk=self.link.pk)

        self.assertEqual(response.status_code, 200, f"Could not update herald")

        pass

    def test_subscribe(self):
        """Test email subscription to herald."""
        # Create fake request
        request = fake_request(url=reverse('herald_subscribe'), data={},user=self.staff_user)
        response = ajax.herald_subscribe(request=request)
        self.assertEqual(response.status_code, 200, f"Could toggle subscription")

        # Toggle subscription
        response = ajax.herald_subscribe(request=request)
        self.assertEqual(response.status_code, 200, f"Could toggle subscription")

        pass
