import logging
import hjson

from django.test import TestCase, RequestFactory
from unittest.mock import patch, MagicMock
from django.urls import reverse
from django.conf import settings

from biostar.engine import auth, const
from biostar.engine import models, views, api

from . import util


logger = logging.getLogger('engine')


class RecipeViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com", is_staff=True)
        self.owner.set_password("test")
        self.factory = RequestFactory()

        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary",
                                           uid="test")
        #Test data
        self.recipe = auth.create_analysis(project=self.project, json_text="{}", template="")
        self.recipe.save()


    @patch('biostar.engine.models.Job.save', MagicMock(name="save"))
    def test_recipe_run(self):
        "Test the recipe run view with POST request"

        data = {"name": "name of the job"}

        url = reverse('recipe_run', kwargs=dict(uid=self.recipe.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        self.recipe.security = models.Analysis.AUTHORIZED
        self.recipe.save()

        response = views.recipe_run(request=request, uid=self.recipe.uid)

        self.process_response(response=response, data=data, save=True, model=models.Job)


    @patch('biostar.engine.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_code_edit(self):
        "Test the recipe preview/save code view with POST request"

        data = {'action': "SAVE", 'template':'#test change', 'json':"{}"}
        url = reverse('recipe_code_edit', kwargs=dict(uid=self.recipe.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_code_edit(request=request, uid=self.recipe.uid)

        self.process_response(response=response, data=data, save=True)


    @patch('biostar.engine.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_edit(self):
        "Test recipe edit with POST request"

        data = { "name": "test", "sticky":True, "summary":"summary", "text":"text" ,
                 "uid":"test"}
        url = reverse('recipe_edit', kwargs=dict(uid=self.recipe.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_edit(request=request, uid=self.recipe.uid)
        self.process_response(response=response, data=data, save=True)

    def test_recipe_copy(self):
        "Test recipe copy interface"

        url = reverse('recipe_copy', kwargs=dict(uid=self.recipe.uid))

        request = util.fake_request(url=url, data={}, user=self.owner)

        response = views.recipe_copy(request=request, uid=self.recipe.uid)

        self.process_response(response=response, data={})

    def test_recipe_paste(self):
        "Test recipe paste interface"

        url = reverse('recipe_paste', kwargs=dict(uid=self.recipe.project.uid))

        request = util.fake_request(url=url, data={}, user=self.owner)

        request.session[settings.CLIPBOARD_NAME] = {const.RECIPE_CLIPBOARD: self.recipe.uid}

        response = views.recipe_paste(request=request, uid=self.recipe.project.uid)

        self.process_response(response=response, data={})

    def Xtest_api(self):
        "Test the recipe api"

        api_list = reverse('recipe_api_list'), api.recipe_api_list, {}
        api_json = reverse('recipe_api_json', kwargs=dict(uid=self.recipe.uid)), api.recipe_json, dict(uid=self.recipe.uid)
        api_template = reverse('recipe_api_template', kwargs=dict(uid=self.recipe.uid)), api.recipe_template, dict(uid=self.recipe.uid)

        for data in [api_list, api_json, api_template]:
            url, view_func, params = data

            request = util.fake_request(url=url, data={'k': settings.API_KEY}, user=self.owner)

            response = view_func(request=request, **params)

            self.assertEqual(response.status_code, 200, f"Could not redirect :\nresponse:{response}")

    def test_recipe_update(self):
        "Test updating recipe through auth"

        changed = auth.create_analysis(project=self.project,
                                       json_text=self.recipe.json_text,
                                       template=self.recipe.template,
                                       uid=self.recipe.uid, update=True)

        self.assertEqual(changed.uid, self.recipe.uid)

    def process_response(self, response, data, model=models.Analysis,save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect :\nresponse:{response}")

        if save:
            self.assertTrue( model.save.called, "save() method not called when editing.")
