import logging
import os
from unittest.mock import patch, MagicMock

from django.test import TestCase, override_settings
from django.urls import reverse
from  django.conf import settings
from biostar.recipes import auth
from biostar.recipes import models, views, forms, api
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
class ProjectViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com",
                                                     is_staff=True)
        self.owner.set_password("tested")

        # Set up project to edit

        self.project = auth.create_project(user=self.owner, name="tested", text="Text", summary="summary",
                                           uid="tested")

    @patch('biostar.recipes.models.Project.save', MagicMock(name="save"))
    def test_create_view(self):
        "Test project create view with POST request"

        request = fake_request(url=reverse('project_create'), data={}, user=self.owner)
        response = views.project_create(request)

        self.process_response(response=response, data={}, save=True)

    @patch('biostar.recipes.models.Project.save', MagicMock(name="save"))
    def test_project_delete(self):

        url = reverse('project_delete', kwargs=dict(uid=self.owner.profile.uid))

        request = fake_request(url=url, method='GET', data={}, user=self.owner)

        response = views.project_delete(request, uid=self.project.uid)

        self.process_response(response, data={})

        return

    def Xtest_api_project_info(self):
        """
        Test replace project info using PUT request
        """
        url = reverse('project_api_info', kwargs=dict(uid=self.project.uid))

        request = fake_request(url=url, method='PUT', data={'k': settings.API_KEY}, user=self.owner)

        response = api.project_info(request, uid=self.project.uid)

        self.process_response(response, data={})

        return

    @patch('biostar.recipes.models.Project.save', MagicMock(name="save"))
    def test_edit_view(self):
        "Test project edit view with POST request"

        # Create fake request
        data = {'name': 'New Name', 'uid': 'tested', "summary": "summary", "rank": 100,
                'text': 'tested', "privacy": models.Project.SHAREABLE}

        url = reverse('project_edit', kwargs=dict(uid=self.project.uid))

        request = fake_request(url=url, data=data, user=self.owner)

        response = views.project_edit(request, uid=self.project.uid)

        self.process_response(response=response, data=data, save=True)

    def test_users_view(self):
        "Test project_users with POST request"

        new_user = models.User.objects.create_user(username="test2", email="test2@l.com")
        new_user.set_password("test2")

        data = {"access": models.Access.WRITE_ACCESS,
                "user_id": new_user.id, "project_uid": self.project.uid}

        url = reverse('project_users', kwargs=dict(uid=self.project.uid))

        request = fake_request(url=url, data=data, user=self.owner)

        response = views.project_users(request, uid=self.project.uid)

        data['uid'] = self.project.uid
        self.assertEqual(response.status_code, 200, f"Error with response code :\nresponse:{response}")

    def test_project_share(self):
        "Test project share by access it using a sharable token."

        url = reverse('project_share', kwargs=dict(token=self.project.sharable_token))

        self.project.privacy = models.Project.SHAREABLE
        self.project.save()

        request = fake_request(url=url, method='GET', data={}, user=self.owner)
        response = views.project_share(request, token=self.project.sharable_token)

        self.process_response(response, data={})

        pass

    def test_project_update(self):
        "Test updating the project"

        changed = auth.create_project(user=self.owner, name="New name",
                                      uid=self.project.uid, update=True)

        self.assertEqual(changed.uid, self.project.uid)

    def process_response(self, response, data, save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to project view after tested :\nresponse:{response}")

        if save:
            self.assertTrue(models.Project.save.called, "save() method not called when editing.")
