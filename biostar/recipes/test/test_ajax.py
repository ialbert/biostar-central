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

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))

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
        """
        Test AJAX function to check and update job status
        """
        data = {'state': models.Job.RUNNING}

        url = reverse('ajax_check_job', kwargs=dict(uid=self.job.uid))

        request = util.fake_request(url=url, data=data, user=self.owner, method='GET')

        json_response = ajax.check_job(request=request, uid=self.job.uid)
        self.process_response(json_response)

    def test_snippet_code(self):
        """
        Test AJAX function used to add snippets to scripts.
        """
        data = {'command': self.snippet.uid, 'template': '# current recipe.'}

        url = reverse('snippet_code')

        request = util.fake_request(url=url, data=data, user=self.owner)

        json_response = ajax.snippet_code(request=request)

        self.process_response(json_response)
        return

    def test_snippet_form(self):
        """
        Test function that renders the snippet create/edit form.
        """
        # Populate form with category data
        category_data = {'is_category': int(True), 'type_name': self.snippet_type.name, 'type_uid': self.snippet_type.uid}

        # Populate form with snippet data.
        snippet_data = {'is_category': int(False), 'snippet_uid': self.snippet.uid,
                        'snippet': self.snippet.command, 'help_text': self.snippet.help_text}

        url = reverse('snippet_form')

        # Test rendering category form
        request = util.fake_request(url=url, data=category_data, user=self.owner)
        json_response = ajax.snippet_form(request=request)
        self.process_response(json_response)

        # Test rendering snippet form
        request = util.fake_request(url=url, data=snippet_data, user=self.owner)
        json_response = ajax.snippet_form(request=request)

        self.process_response(json_response)

    def test_create_snippet_type(self):
        """
        Test function that creates snippet types ( categories )
        """

        path = os.path.join(CURRENT_DIR, "data", "image.png")
        image_stream = open(path, "rb")

        data = {'name': self.snippet.uid, 'image': image_stream}

        url = reverse('create_snippet_type')

        request = util.fake_request(url=url, data=data, user=self.owner)

        json_response = ajax.create_snippet_type(request=request)

        self.process_response(json_response)
        image_stream.close()

    def test_create_snippet(self):
        """
        Test function that creates snippets
        """
        create_data = {'snippet': 'ls -l', 'help_text': ' Help text',
                       'type_uid': self.snippet_type.uid}

        url = reverse('create_snippet')

        # Test creating a snippet
        request = util.fake_request(url=url, data=create_data, user=self.owner)
        json_response = ajax.create_snippet(request=request)

        self.process_response(json_response)

    def test_edit_snippet(self):
        """
        Test AJAX function that edits snippets.
        """

        edit_data = {'snippet_uid': self.snippet.uid, 'snippet': 'ls -l', 'help_text': ' Help text',
                     'type_uid': self.snippet_type.uid}

        url = reverse('create_snippet')

        # Test creating a snippet
        request = util.fake_request(url=url, data=edit_data, user=self.owner)
        json_response = ajax.create_snippet(request=request)

        self.process_response(json_response)

    def test_preview_template(self):
        """
        Test AJAX function used to preview recipe scripts
        """

        data = {'template': '# recipe code', 'json_text': '{}', 'name': self.recipe.name,
                'uid': self.recipe.uid}

        url = reverse('preview_template')
        request = util.fake_request(url=url, data=data, user=self.owner)
        json_response = ajax.preview_template(request=request)

        self.process_response(json_response)

    def test_preview_json(self):
        """
        Test AJAX function used to preview recipe json
        """

        data = {'name': self.recipe.name, 'project_uid':self.recipe.project.uid,
                'json_text': self.recipe.json_text}

        url = reverse('preview_json')
        request = util.fake_request(url=url, data=data, user=self.owner)
        json_response = ajax.preview_json(request=request)

        self.process_response(json_response)

    def test_add_recipe_fields(self):
        """
        Test AJAX function used to add fields to the recipe JSON
        """

        display_types = ['radio', 'data', 'integer', 'textbox', 'float', 'checkbox',
                         'dropdown']

        json_text = '{}'
        url = reverse('add_recipe_fields')

        for dtype in display_types:
            data = {'json_text': json_text, 'display_types': dtype}

            request = util.fake_request(url=url, data=data, user=self.owner)
            json_response = ajax.add_to_interface(request=request)
            self.process_response(json_response)

    def process_response(self, response):
        "Check the response on POST request is redirected"

        response_data = response.content
        response_data = json.loads(response_data)

        self.assertEqual(response_data['status'], 'success', f'Error :{response_data["msg"]}')




