
import logging, os
from django.test import TestCase
from unittest.mock import patch, MagicMock
from django.urls import reverse
from django.conf import settings

from biostar.engine import auth
from biostar.engine import models, views, forms

from . import util

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'engine', 'test'))

logger = logging.getLogger('engine')



class ProjectViewTest(TestCase):


    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com")
        self.owner.set_password("test")

        # Set up project to edit

        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary",
                                           uid="testing")

    @patch('biostar.engine.models.Project.save', MagicMock(name="save"))
    def test_create_view(self):
        "Test project create view with POST request"

        path = os.path.join(TEST_ROOT, "data", "image.png")
        image_stream = open(path, "rb")

        # Create fake request
        data = {'name': 'My project', 'uid': 'example', "summary":"summary",
                'text': 'testing', "privacy": models.Project.PRIVATE, "image":image_stream}

        request = util.fake_request(url=reverse('project_create'), data=data, user=self.owner)
        response = views.project_create(request)

        image_stream.close()

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

        self.process_response(response=response, data=data, save=True)


    def test_users_view(self):
        "Test project_users with POST request"

        new_user = models.User.objects.create_user(username="test2", email="test2@l.com")
        new_user.set_password("test2")

        data = {"access":models.Access.WRITE_ACCESS,
                "user_id":new_user.id, "project_uid":self.project.uid}

        url = reverse('project_users', kwargs=dict(uid=self.project.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.project_users(request, uid=self.project.uid)

        data['uid'] = self.project.uid
        self.process_response(response=response, data=data)


    def test_project_update(self):
        "Test updating the project"

        changed = auth.create_project(user=self.owner, name="New name",
                                      uid=self.project.uid, update=True)

        self.assertEqual(changed.uid, self.project.uid)

    def test_access_forms(self):
        " Test generating a list of forms for a given user list"

        new_user = models.User.objects.create_user(username="test2", email="test2@l.com")
        new_user.set_password("test2")

        users = [self.owner, new_user ]

        user_forms = forms.access_forms(users, project=self.project)

        # Error generating users access forms ( forms.access_forms).
        self.assertTrue(len(users) ==len(user_forms))

    def process_response(self, response, data, save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to project view after testing :\nresponse:{response}")

        if save:
            self.assertTrue( models.Project.save.called, "save() method not called when editing.")

