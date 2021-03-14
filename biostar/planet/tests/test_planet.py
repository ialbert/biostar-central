import logging
import json
import os
from django.test import TestCase
from django.urls import reverse
from django.test import TestCase, override_settings

from biostar.accounts.models import User, Profile

from biostar.planet import models, views, auth
from biostar.utils.helpers import fake_request
from biostar.forum.util import get_uuid

__MODULE_DIR = os.path.dirname(auth.__file__)
TEST_ROOT = os.path.join(__MODULE_DIR, 'tests')

logger = logging.getLogger('engine')

PLANET_DIR = os.path.abspath(os.path.join(TEST_ROOT, "feeds"))


@override_settings(PLANET_DIR=PLANET_DIR, INIT_PLANET=False)
class PlanetTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"test", email="tested@tested.com", password="tested",
                                         is_superuser=True)
        pass

    def test_blog_create(self):
        """
        Test Planet blog creation.
        """
        fname = os.path.join(PLANET_DIR, 'test-feeds.txt')
        auth.add_blogs(add_fname=fname)

        feed = open(fname, 'r').readline()

        # Check to see if the blog post has been created
        blog_query = models.Blog.objects.filter(feed=feed)

        blog = blog_query.first()
        print(blog.fname)
        blog.download()
        self.assertTrue(blog_query.exists(), "Error creating blog in database")
        self.assertTrue(blog.parse() is not None, "Error parsing blog")
