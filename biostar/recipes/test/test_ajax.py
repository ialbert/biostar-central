import logging
import os
import json
from unittest.mock import patch, MagicMock

from django.conf import settings
from django.test import TestCase, override_settings
from django.urls import reverse

from biostar.recipes import models, auth, ajax
from biostar.utils.helpers import fake_request, get_uuid

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'tested'))
TOC_ROOT = os.path.join(TEST_ROOT, 'toc')

# Ensure that the table of directory exists.
os.makedirs(TOC_ROOT, exist_ok=True)

__MODULE_DIR = os.path.dirname(auth.__file__)
TEST_DIR = os.path.join(__MODULE_DIR, 'test')

IMPORT_ROOT_DIR = os.path.join(TEST_DIR, 'data')
logger = logging.getLogger('engine')


@override_settings(MEDIA_ROOT=TEST_ROOT, TOC_ROOT=TOC_ROOT, IMPORT_ROOT_DIR=IMPORT_ROOT_DIR)
class AjaxTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com")
        self.owner.set_password("tested")

        # Set up generic owner
        self.trusted_owner = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@2.com",
                                                             is_superuser=True)
        self.owner.set_password("tested")

        self.project = auth.create_project(user=self.owner, name="tested", text="Text", summary="summary",
                                           uid="tested")
        self.project2 = auth.create_project(user=self.trusted_owner, name="tested", text="Text", summary="summary",
                                            uid="tested")

        self.recipe = auth.create_analysis(project=self.project, json_text="", template="",
                                           security=models.Analysis.AUTHORIZED)

        self.snippet_type = models.SnippetType.objects.create(name='Snippet type', owner=self.owner)

        self.snippet = models.Snippet.objects.create(command='ls -l', type=self.snippet_type,
                                                     help_text='List files in directory',
                                                     owner=self.owner)

        self.job = auth.create_job(analysis=self.recipe, user=self.owner)
        self.job.save()

    def test_check_job(self):
        """
        Test AJAX function to check and update job status
        """
        data = {'state': models.Job.RUNNING}

        url = reverse('ajax_check_job', kwargs=dict(uid=self.job.uid))

        request = fake_request(url=url, data=data, user=self.owner, method='GET')

        json_response = ajax.check_job(request=request, uid=self.job.uid)
        self.process_response(json_response)

    def test_copy_file(self):
        """
        Test AJAX function used to copy file
        """
        data = {'path': os.path.join(IMPORT_ROOT_DIR, "plain-text.txt")}

        url = reverse('copy_file')

        request = fake_request(url=url, data=data, user=self.trusted_owner)

        json_response = ajax.copy_file(request=request)

        self.process_response(json_response)

    def test_toggle_delete(self):
        """
        Test AJAX function used to toggle delete on objects
        """
        data = {'uid': self.job.uid, "type": 'job'}

        url = reverse('toggle_delete')

        request = fake_request(url=url, data=data, user=self.owner)

        json_response = ajax.toggle_delete(request=request)

        self.process_response(json_response)

    def test_manage_access(self):
        """
        Test AJAX function used to manage user access
        """
        user2 = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com")
        data = {'user_id': user2.id, "project_uid": self.project.uid, "access": 'write'}

        url = reverse('toggle_delete')

        request = fake_request(url=url, data=data, user=self.owner)

        json_response = ajax.manage_access(request)

        self.process_response(json_response)

    def test_copy_object(self):
        """
        Test AJAX function used to copy objects
        """

        user2 = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com")
        data = {'user_id': user2.id, "project_uid": self.project.uid, "access": 'write'}

        url = reverse('toggle_delete')

        request = fake_request(url=url, data=data, user=self.owner)

        json_response = ajax.manage_access(request)

        self.process_response(json_response)

    def test_preview_json(self):
        """
        Test AJAX function used to preview recipe json
        """

        data = {'recipe':self.recipe.id, 'toml': "[foo]\nparam=2"}

        url = reverse('preview_json')
        request = fake_request(url=url, data=data, user=self.owner)
        json_response = ajax.preview_json(request=request)

        self.process_response(json_response)

    def process_response(self, response):
        "Check the response on POST request is redirected"

        response_data = response.content
        response_data = json.loads(response_data)

        self.assertEqual(response_data['status'], 'success', f'Error :{response_data["msg"]}')




