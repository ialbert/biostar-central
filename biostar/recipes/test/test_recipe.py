import logging
import os
from unittest.mock import patch, MagicMock

from django.conf import settings
from django.test import TestCase, RequestFactory, override_settings
from django.urls import reverse
from django.test import TestCase, override_settings
#from biostar.accounts.models import Use

from biostar.recipes import auth, const
from biostar.recipes import models, views, api
from biostar.utils.helpers import fake_request, get_uuid

logger = logging.getLogger('engine')

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'tested'))
TOC_ROOT = os.path.join(TEST_ROOT, 'toc')
__CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))

# Ensure that the table of directory exists.
os.makedirs(TOC_ROOT, exist_ok=True)

@override_settings(MEDIA_ROOT=TEST_ROOT, TOC_ROOT=TOC_ROOT)
class RecipeRunTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com",
                                                     is_staff=True, is_superuser=True)
        self.owner.set_password("tested")
        self.factory = RequestFactory()

        self.project = auth.create_project(user=self.owner, name="tested", text="Text", summary="summary",
                                           uid="tested")
        # Test data
        self.recipe = auth.create_analysis(project=self.project, json_text="", template="#test template")
        self.recipe.save()

    def test_authorize_run(self):
        """Test to see if function that authorizes runs works correctly."""

        user1 = self.owner
        recipe = self.recipe

        # Current user can run the recipe

        self.assertTrue(auth.authorize_run(user1, recipe), "Authorized users can not run recipes.")

        user2 = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com")

        self.assertFalse(auth.authorize_run(user2, recipe), "Unauthorized users can run recipes.")
        return


@override_settings(MEDIA_ROOT=TEST_ROOT)
class RecipeViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com",
                                                     is_staff=True, is_superuser=True)
        self.owner.set_password("tested")
        self.factory = RequestFactory()

        self.project = auth.create_project(user=self.owner, name="tested", text="Text", summary="summary",
                                           uid="tested")
        # Test data
        self.recipe = auth.create_analysis(project=self.project, json_text="", template="#test template")
        self.recipe.save()

    @patch('biostar.recipes.models.Job.save', MagicMock(name="save"))
    def test_recipe_run(self):
        "Test the recipe run view with POST request"

        data = {"name": "name of the job"}

        url = reverse('recipe_run', kwargs=dict(uid=self.recipe.uid))

        request = fake_request(url=url, data=data, user=self.owner)

        self.recipe.security = models.Analysis.AUTHORIZED
        self.recipe.save()

        response = views.recipe_run(request=request, uid=self.recipe.uid)

        self.process_response(response=response, data=data, save=True, model=models.Job)

    @patch('biostar.recipes.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_create(self):
        "Test recipe create with POST request"
        data = {"name": "tested", "summary": "summary", "text": "text", "rank": 100,
                "uid": "tested", 'json_text':'', 'template':'# Code here'}
        url = reverse('recipe_create', kwargs=dict(uid=self.project.uid))

        request = fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_create(request=request, uid=self.project.uid)

        self.process_response(response=response, data=data, save=True)

    @patch('biostar.recipes.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_edit(self):
        "Test recipe edit with POST request"
        from biostar.recipes import ajax

        data = {"name": "tested", "text": "text", "rank": 100,
                "uid": "tested", 'json_text':'', 'template':'# Code here'}
        url = reverse('ajax_recipe_edit', kwargs=dict(id=f"{self.recipe.id}"))

        request = fake_request(url=url, data=data, user=self.owner)

        response = ajax.ajax_edit(request=request, id=self.recipe.id)

        #self.process_response(response=response, data=data, save=True)

    def test_recipe_code_download(self):
        "Test recipe code download "

        url = reverse("recipe_download", kwargs=dict(uid=self.recipe.uid))
        request = fake_request(url=url, data={}, user=self.owner)
        response = views.recipe_code_download(request=request, uid=self.recipe.uid)

        self.assertTrue(response.content.decode() == self.recipe.template,
                        f"Error downloading code. Expected: {self.recipe.template} "
                        f"received: {response.content.decode()}")

    def test_recipe_delete(self):
        "Test reset delete"

        url = reverse('recipe_delete', kwargs=dict(uid=self.recipe.uid))

        request = fake_request(url=url, data={}, user=self.owner)

        response = views.recipe_delete(request=request, uid=self.recipe.uid)

        self.process_response(response=response, data={})

    def Xtest_api(self):
            "Test the recipe api"

            api_list = reverse('api_list'), api.recipe_api_list, {}
            api_json = reverse('recipe_api_json', kwargs=dict(uid=self.recipe.uid)), api.recipe_json, dict(
                uid=self.recipe.uid)
            api_template = reverse('recipe_api_template', kwargs=dict(uid=self.recipe.uid)), api.recipe_template, dict(
                uid=self.recipe.uid)

            for data in [api_list, api_json, api_template]:
                url, view_func, params = data

                request = fake_request(url=url, data={'k': settings.API_KEY}, user=self.owner)

                response = view_func(request=request, **params)

                self.assertEqual(response.status_code, 200, f"Could not redirect :\nresponse:{response}")

    def test_recipe_update(self):
        "Test updating recipe through auth"

        changed = auth.create_analysis(project=self.project,
                                       json_text=self.recipe.json_text,
                                       template=self.recipe.template,
                                       uid=self.recipe.uid, update=True)

        self.assertEqual(changed.uid, self.recipe.uid)

    def process_response(self, response, data, model=models.Analysis, save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect :\nresponse:{response}")

        if save:
            self.assertTrue(model.save.called, "save() method not called when editing.")
