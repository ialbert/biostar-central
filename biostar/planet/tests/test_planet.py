import logging
import json
from django.test import TestCase
from django.urls import reverse

from biostar.accounts.models import User, Profile

from biostar.planet import models, views, auth
from biostar.utils.helpers import fake_request
from biostar.forum.util import get_uuid

logger = logging.getLogger('engine')


class PlanetTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"test", email="tested@tested.com", password="tested",
                                         is_superuser=True)
        pass

    def Xtest_blog_create(self):

        for stype in ["unfollow", "messages", "email", "all", "default"]:

            data = {"sub_type": stype, "root_uid": self.post.uid}
            request = fake_request(url=reverse('vote'), data=data, user=self.owner)
            response = ajax.ajax_subs(request)
            self.assertEqual(response.status_code, 200, f"Could not preform subscription action:{stype}.")
