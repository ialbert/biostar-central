import logging
import os
import shutil
import datetime
from django.core import management
from django.urls import reverse
from django.test import TestCase, override_settings
from django.conf import settings
from biostar.forum import models, api
from biostar.utils.helpers import fake_request
from biostar.accounts.models import User

logger = logging.getLogger('engine')

TEST_DATABASE_NAME = f"test_{settings.DATABASE_NAME}"
TEST_DEBUG = True

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'test'))
TEST_INDEX_DIR = TEST_ROOT
TEST_INDEX_NAME = "index"



class PostAPITest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"test", email="tested@tested.com", password="tested")
        self.staff_user = User.objects.create(username=f"test2", is_superuser=True, is_staff=True,
                                              email="tested@staff.com", password="tested")

        # Create an existing tested post
        self.post = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                     type=models.Post.QUESTION)
        self.owner.save()
        pass

    def test_api(self):
        """Test post creation with POST request"""

        url = reverse("api_stats_on_day", kwargs=dict(day=0))

        request = fake_request(url=url, data={}, user=self.owner)

        response = api.daily_stats_on_day(request=request, day=0)

        print(response, response.content)
        self.assertEqual(response.status_code, 404)

        date = datetime.datetime.now() - datetime.timedelta(days=1)
        stats = api.compute_stats(date)
        print(date, stats)

        url = reverse("api_stats_on_date", kwargs=dict(year=2020, month=10, day=1))

        request = fake_request(url=url, data={}, user=self.owner)

        response = api.daily_stats_on_date(request=request, year=2020, month=10, day=1)
        self.assertEqual(response.status_code, 200)
        #self.process_response(response=response)


