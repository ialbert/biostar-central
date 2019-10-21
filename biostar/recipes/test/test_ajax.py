import logging
import os
import json
from unittest.mock import patch, MagicMock

from django.conf import settings
from django.test import TestCase, override_settings
from django.urls import reverse

from biostar.recipes import models, auth, ajax
from . import util

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'tested'))

__CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))

logger = logging.getLogger('engine')


@override_settings(MEDIA_ROOT=TEST_ROOT)
class AjaxTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username=f"tested{util.get_uuid(10)}", email="tested@l.com")
        self.owner.set_password("tested")

        self.project = auth.create_project(user=self.owner, name="tested", text="Text", summary="summary",
                                           uid="tested")

        self.recipe = auth.create_analysis(project=self.project, json_text="{}", template="",
                                           security=models.Analysis.AUTHORIZED)

        self.snippet_type = models.SnippetType.objects.create(name='Snippet type', owner=self.owner)

        self.snippet = models.Snippet.objects.create(command='ls -l', type=self.snippet_type,
                                                     help_text='List files in directory',
                                                     owner=self.owner)

        self.job = auth.create_job(analysis=self.recipe, user=self.owner)
        self.job.save()

    def test_check_job(self):

        data = {'state': models.Job.RUNNING}

        url = reverse('ajax_check_job', kwargs=dict(uid=self.job.uid))

        request = util.fake_request(url=url, data=data, user=self.owner, method='GET')

        json_response = ajax.check_job(request=request, uid=self.job.uid)

        response_data = json_response.content
        response_data = json.loads(response_data)

        self.assertEqual(response_data['status'], 'success', f'Error checking job:{response_data["msg"]}')

    def test_snippet_code(self):

        data = {'command': self.snippet.uid, 'template': '# current recipe.'}

        url = reverse('snippet_code')

        request = util.fake_request(url=url, data=data, user=self.owner)

        json_response = ajax.snippet_code(request=request)

        response_data = json_response.content
        response_data = json.loads(response_data)

        print(response_data)
        1/0
        self.assertEqual(response_data['status'], 'success', f'Error adding snippet:{response_data["msg"]}')

        return

    def test_snippet_form(self):

        return


    def test_create_snippet_type(self):
        return


