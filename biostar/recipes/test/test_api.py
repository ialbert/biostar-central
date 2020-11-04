import logging
import  json
import os
from unittest.mock import patch, MagicMock

from django.test import TestCase, override_settings
from django.urls import reverse
from  django.conf import settings
from biostar.recipes import api, auth, models
from biostar.utils.helpers import fake_request, get_uuid

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'tested'))
TOC_ROOT = os.path.join(TEST_ROOT, 'toc')
__MODULE_DIR = os.path.dirname(auth.__file__)
TEST_DIR = os.path.join(__MODULE_DIR, 'test')

__CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
logger = logging.getLogger('engine')

# Ensure that the table of directory exists.
os.makedirs(TOC_ROOT, exist_ok=True)


@override_settings(MEDIA_ROOT=TEST_ROOT, TOC_ROOT=TOC_ROOT)
class APITest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com",
                                                     is_staff=True)
        self.owner.set_password("tested")

        # Set up project to edit

        self.project = auth.create_project(user=self.owner, name="tested", text="Text", summary="summary",
                                           uid="tested")
        # Set up generic data for editing
        self.data = auth.create_data(project=self.project, path=__file__, name="tested")

    @patch('biostar.recipes.models.Project.save', MagicMock(name="save"))
    def test_project_api(self):
        "Test project api"

        # Edit project
        newuid = self.project.uid
        view = reverse('project_api')
        fname = os.path.join(TEST_DIR, "data", "demo.json")
        stream = open(fname, 'r')

        data = dict(token=self.owner.profile.token, uid=newuid, data=stream)
        request = fake_request(url=view, data=data, user=self.owner, method="POST")

        response = api.project_api(request=request)
        data = json.loads(response.content.decode())
        self.assertEqual(response.status_code, 200)

    @patch('biostar.recipes.models.Project.save', MagicMock(name="save"))
    def test_recipe_api(self):

        # Edit project
        rec = self.project.analysis_set.first()
        newuid = rec.uid
        view = reverse('recipe_api')
        fname = os.path.join(TEST_DIR, "data", "demo.json")
        stream = open(fname, 'r')

        data = dict(token=self.owner.profile.token,
                    uid=newuid, pid=self.project.uid,
                    data=stream)
        request = fake_request(url=view, data=data, user=self.owner, method="POST")

        response = api.recipe_api(request=request)
        #data = json.loads(response.content.decode())
        self.assertEqual(response.status_code, 200)

        return

    def test_data_api(self):
        """
        Test data api
        """
        url = reverse('data_api')
        data = dict(uid=self.data.uid, token=self.owner.profile.token)
        request = fake_request(url=url, data=data, user=self.owner, method="GET")
        response = api.data_api(request=request)
        self.assertEqual(response.status_code, 200)

        return

    def process_response(self, response):
        "Check the response on POST request is redirected"

        response_data = response.content
        response_data = json.loads(response_data)

        self.assertEqual(response_data['status'], 'success', f'Error :{response_data["msg"]}')



