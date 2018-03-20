import logging
import hjson

from django.test import TestCase, RequestFactory
from unittest.mock import patch, MagicMock
from django.urls import reverse

from biostar.engine import auth, const
from biostar.engine import models, views

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
                                           uid="testing")
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
    def test_recipe_code(self):
        "Test the recipe preview/save code view with POST request"

        data = {'action': "SAVE", 'template':'#test change', 'json':"{}"}
        url = reverse('recipe_code', kwargs=dict(uid=self.recipe.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_code(request=request, uid=self.recipe.uid)

        self.process_response(response=response, data=data, save=True)

    @patch('biostar.engine.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_create(self):
        "Test recipe create with POST request"

        data = { "name": "test", "sticky":True, "summary":"summary", "text":"text" }
        url = reverse('recipe_create', kwargs=dict(uid=self.project.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_create(request=request, uid=self.project.uid)

        self.process_response(response=response, data=data, save=True)

    @patch('biostar.engine.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_edit(self):
        "Test recipe edit with POST request"

        data = { "name": "test", "sticky":True, "summary":"summary", "text":"text" ,
                 "uid":"testing"}
        url = reverse('recipe_edit', kwargs=dict(uid=self.recipe.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_edit(request=request, uid=self.recipe.uid)
        self.process_response(response=response, data=data, save=True)


    def test_recipe_mod_page(self):
        "test moderator view works"

        url = reverse('recipe_mod')
        request = util.fake_request(url=url, data={}, user=self.owner)

        response = views.recipe_mod(request=request)

        self.assertEqual(response.status_code, 200, "Can not load moderate page.")


    @patch('biostar.engine.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_diff(self):
        "Test the recipe diff generation."

        self.recipe.last_valid = "Test"
        self.recipe.save()
        url = reverse('recipe_diff', kwargs=dict(uid=self.recipe.uid))

        approve_request = util.fake_request(url=url, data={'action':"APPROVE"}, user=self.owner)
        revert_request = util.fake_request(url=url, data={'action':"REVERT"}, user=self.owner)

        for action in (approve_request, revert_request):

            response = views.recipe_diff(request=action, uid=self.recipe.uid)

            self.process_response(response=response, data={'action':action}, save=True)


    def process_response(self, response, data, model=models.Analysis,save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect :\nresponse:{response}")

        if save:
            self.assertTrue( model.save.called, "save() method not called when editing.")