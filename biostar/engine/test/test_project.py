
import logging, os
from django.test import TestCase
from unittest.mock import patch, MagicMock
from django.urls import reverse
from django.contrib.messages.storage import fallback

from biostar.engine import auth
from biostar.engine import models, views, forms

from . import util

logger = logging.getLogger('engine')

class ProjectViewTest(TestCase):


    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com")
        self.owner.set_password("test")

        # Set up project to edit
        pre = models.Project.objects.count()
        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary",
                                           uid="testing")

        self.assertTrue(models.Project.objects.count() == (pre + 1), "Error creating project in database")

    @patch('biostar.engine.models.Project.save', MagicMock(name="save"))
    def test_create_view(self):
        "Test project create view with POST request"

        # Create fake request
        data = {'name': 'My project', 'uid': 'example', "summary":"summary",
                'text': 'testing', "privacy": models.Project.PRIVATE}

        request = util.fake_request(url=reverse('project_create'), data=data, user=self.owner)

        response = views.project_create(request)

        util.remove_test_folders(self.project.get_project_dir())

        self.process_response(response=response, data=data, save=True)


    @patch('biostar.engine.models.Project.save', MagicMock(name="save"))
    def test_edit_view(self):
        "Test project edit view with POST request"

        # Create fake request
        data = {'name': 'New Name', 'uid': 'testing', "summary":"summary",
                'text': 'testing', "privacy": models.Project.SHAREABLE}

        url = reverse('project_edit', kwargs=dict(uid=self.project.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.project_edit(request, uid=self.project.uid)

        util.remove_test_folders(self.project.get_project_dir())

        self.process_response(response=response, data=data, save=True)


    def test_users_view(self):
        "Test project_users with POST request"

        new_user = models.User.objects.create_user(username="test2", email="test2@l.com")
        new_user.set_password("test2")

        data = {"access":models.Access.ADMIN_ACCESS,
                "user_id":new_user.id, "project_id":self.project.id}

        url = reverse('project_users', kwargs=dict(uid=self.project.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.project_users(request, uid=self.project.uid)

        util.remove_test_folders(self.project.get_project_dir())

        data['uid'] = self.project.uid
        self.process_response(response=response, data=data)


    def test_access_forms(self):
        " Test generating a list of forms for a given user list"

        new_user = models.User.objects.create_user(username="test2", email="test2@l.com")
        new_user.set_password("test2")

        users = [self.owner, new_user ]

        user_forms = forms.access_forms(users, project=self.project)

        # Still need to call this since setup creates a dir
        util.remove_test_folders(self.project.get_project_dir())

        self.assertTrue(len(users) ==len(user_forms), "Error generating users access forms ( forms.access_forms) ")


    def process_response(self, response, data, save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to project view after testing :\nresponse:{response}")

        self.assertTrue(data['uid'] in response.url,
                        "Was not redirected to the correct project.")
        if save:
            self.assertTrue( models.Project.save.called,
                         "project.save() method not called when editing.")

