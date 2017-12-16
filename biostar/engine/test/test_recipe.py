import logging
from django.test import TestCase, RequestFactory
from unittest.mock import patch, MagicMock
from django.urls import reverse

from biostar.engine import auth
from biostar.engine import models, views

from . import util


logger = logging.getLogger('engine')



class RecipeViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com")
        self.owner.set_password("test")
        self.factory = RequestFactory()

        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary",
                                           uid="testing")

        pre = models.Analysis.objects.count()
        self.recipe = auth.create_analysis(project=self.project, json_text="{}", template="")
        self.recipe.save()
        self.assertTrue(models.Analysis.objects.count() == (pre + 1), "Error creating Analysis in database")


    def test_recipe_view(self):
        "Test the recipe copy view with POST request"

        data = {"project":0}

        url = reverse('recipe_view', kwargs=dict(id=self.recipe.id))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_view(request=request, id=self.recipe.id)

        self.assertEqual(response.status_code, 200,
                         f"Could not load to after copying:\nresponse:{response}")

    @patch('biostar.engine.models.Job.save', MagicMock(name="save"))
    def test_recipe_run(self):
        "Test the recipe run view with POST request"

        data = {"name": "name of the job"}

        url = reverse('analysis_run', kwargs=dict(id=self.recipe.id))

        request = util.fake_request(url=url, data=data, user=self.owner)

        self.recipe.security = models.Analysis.AUTHORIZED
        self.recipe.save()

        response = views.recipe_run(request=request, id=self.recipe.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to after running recipe:\nresponse:{response}")

        self.assertTrue("job/list/" in response.url,
                        f"Could not redirect to job list after running recipe:\nresponse:{response}")

        self.assertTrue( models.Job.save.called, "job.save() method not called when running analysis.")


    @patch('biostar.engine.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_code(self):
        "Test the recipe preview/save code view with POST request"

        data = {'action': "SAVE", 'template':'#test change', 'json':'{}'}
        url = reverse('recipe_code', kwargs=dict(id=self.recipe.id))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_code(request=request, id=self.recipe.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to after editing code:\nresponse:{response}")

        self.assertTrue(self.recipe.url() == response.url,
                        f"Could not redirect to correct page: {self.recipe.url()} != {response.url}")

        self.assertTrue( models.Analysis.save.called, "analysis.save() method not called when saving in views.")


    @patch('biostar.engine.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_create(self):
        "Test recipe create with POST request"

        data = { "name": "test", "sticky":True, "summary":"summary", "text":"text" }
        url = reverse('recipe_create', kwargs=dict(uid=self.project.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_create(request=request, uid=self.project.uid)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect after creating recipe:\nresponse:{response}")

        self.assertTrue(f"/recipe/list/{self.project.uid}/" == response.url,
                        f"Could not redirect to correct page: 'recipe/view' != {response.url}")

        self.assertTrue( models.Analysis.save.called, "analysis.save() method not called when creating.")


    @patch('biostar.engine.models.Analysis.save', MagicMock(name="save"))
    def test_recipe_edit(self):
        "Test recipe edit with POST request"

        data = { "name": "test", "sticky":True, "summary":"summary", "text":"text" }
        url = reverse('recipe_edit', kwargs=dict(id=self.recipe.id))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.recipe_edit(request=request, id=self.recipe.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect after creating recipe:\nresponse:{response}")

        self.assertTrue( models.Analysis.save.called, "analysis.save() method not called when creating.")
