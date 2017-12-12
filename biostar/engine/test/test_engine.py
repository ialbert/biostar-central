import hjson, logging, os
from django import forms
from django.core import management
from django.test import TestCase, RequestFactory
from unittest.mock import patch, MagicMock
from django.urls import reverse
from django.contrib.messages.storage import fallback

from biostar.engine import auth, factory
from biostar.engine import models, views

logger = logging.getLogger('engine')



class ProjectViewTest(TestCase):


    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com")
        self.owner.set_password("test")
        self.factory = RequestFactory()

        # Set up project to edit
        pre = models.Project.objects.count()
        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary")

        self.assertTrue(models.Project.objects.count() == (pre + 1), "Error creating project in database")


    @patch('biostar.engine.models.Project.save', MagicMock(name="save"))
    def test_create_view(self):
        "Test project create view with POST request"

        # Create fake request
        data = {'name': 'My project', 'uid': 'example', "summary":"summary",
                'text': 'testing', "privacy": models.Project.PRIVATE}

        request = self.factory.post(reverse('project_create'), data)
        request.session = {}
        messages = fallback.FallbackStorage(request=request)
        request._messages = messages
        request.user = self.owner

        response = views.project_create(request)

        self.process_response(response=response, data=data, save=True)


    @patch('biostar.engine.models.Project.save', MagicMock(name="save"))
    def test_edit_view(self):
        "Test project edit view with POST request"

        # Create fake request
        data = {'name': 'New Name', 'uid': 'new', "summary":"summary",
                'text': 'testing', "privacy": models.Project.SHAREABLE}

        request = self.factory.post(reverse('project_edit', kwargs=dict(uid=self.project.uid)),
                                    data)
        request.session = {}
        messages = fallback.FallbackStorage(request=request)
        request._messages = messages
        request.user = self.owner

        response = views.project_edit(request, uid=self.project.uid)

        self.process_response(response=response, data=data, save=True)


    def test_users_view(self):
        "Test project_users with POST request"

        new_user = models.User.objects.create_user(username="test2", email="test2@l.com")
        new_user.set_password("test2")

        data = {"access":models.Access.ADMIN_ACCESS,
                "user_id":new_user.id, "project_id":self.project.id}
        request = self.factory.post(reverse('project_users', kwargs=dict(uid=self.project.uid)),
                                    data)
        request.session = {}
        messages = fallback.FallbackStorage(request=request)
        request._messages = messages
        request.user = self.owner
        response = views.project_users(request, uid=self.project.uid)

        data['uid'] = self.project.uid
        self.process_response(response=response, data=data)


    def process_response(self, response, data, save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to project view after testing :\nresponse:{response}")

        self.assertTrue(data['uid'] in response.url,
                        "Was not redirected to the correct project.")
        if save:
            self.assertTrue( models.Project.save.called,
                         "project.save() method not called when editing.")


class DataViewTest(TestCase):


    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com")
        self.owner.set_password("test")
        self.factory = RequestFactory()

        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary")

        # Set up generic data for editing
        pre = models.Data.objects.count()
        self.data = auth.create_data(project=self.project, path=__file__)
        self.assertTrue(models.Data.objects.count() == (pre + 1), "Error creating Data in database")

    def test_data_copy_view(self):

        "Test Data copy (create a new project and copy) in views with POST request"

        # 0 is the option picked whe creating and copying
        data = {"project":0}

        request = self.factory.post(reverse('data_view', kwargs=dict(id=self.data.id)),
                                    data)
        request.session = {}
        messages = fallback.FallbackStorage(request=request)
        request._messages = messages
        request.user = self.owner
        response = views.data_view(request=request, id=self.data.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to data view after copying Data:\nresponse:{response}")

    @patch('biostar.engine.models.Data.save', MagicMock(name="save"))
    def test_data_edit(self):
        "Test Data edit view with POST request"

        data = {'name':"new_data", 'summary':"summary", 'text':"testing",
                'sticky':True}

        request = self.factory.post(reverse('data_edit', kwargs=dict(id=self.data.id)),
                                    data)
        request.session = {}
        messages = fallback.FallbackStorage(request=request)
        request._messages = messages
        request.user = self.owner

        response = views.data_edit(request=request, id=self.data.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to data view after copying Data:\nresponse:{response}")

        self.assertTrue( models.Data.save.called, "data.save() method not called when editing.")


    @patch('biostar.engine.models.Data.save', MagicMock(name="save"))
    def test_data_upload(self):
        "Test Data upload POST request"

        data = {'file':__file__, 'summary':'summary', "text":"testing", "sticky":True}

        request = self.factory.post(reverse('data_upload', kwargs=dict(uid=self.project.uid)),
                                    data)
        request.session = {}
        messages = fallback.FallbackStorage(request=request)
        request._messages = messages
        request.user = self.owner

        #TODO: this is not showing up in coverage wtffff
        response = views.data_upload(request=request, uid=self.project.uid)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to after uploading:\nresponse:{response}")
        self.assertTrue( "data/list/" in response.url,
                         f"Could not redirect to data list after uploading:\nresponse:{response}")



class RecipeViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com")
        self.owner.set_password("test")
        self.factory = RequestFactory()

        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary")

        pre = models.Analysis.objects.count()
        self.recipe = auth.create_analysis(project=self.project, json_text="{}", template="")
        self.assertTrue(models.Analysis.objects.count() == (pre + 1), "Error creating Analysis in database")


    def test_recipe_view(self):
        "Test the recipe copy view with POST request"

        data = {"project":0}

        request = self.factory.post(reverse('recipe_view', kwargs=dict(id=self.recipe.id)),
                                    data)
        request.session = {}
        messages = fallback.FallbackStorage(request=request)
        request._messages = messages
        request.user = self.owner

        response = views.recipe_view(request=request, id=self.recipe.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to after copying:\nresponse:{response}")

    def test_recipe_run(self):
        "Test the recipe run view with POST request"

        data = {"name": "name of the job"}

        request = self.factory.post(reverse('analysis_run', kwargs=dict(id=self.recipe.id)),
                                    data)
        request.session = {}
        messages = fallback.FallbackStorage(request=request)
        request._messages = messages
        request.user = self.owner

        self.recipe.security = models.Analysis.AUTHORIZED
        self.recipe.save()

        response = views.recipe_run(request=request, id=self.recipe.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to after running recipe:\nresponse:{response}")

        self.assertTrue("job/list/" in response.url,
                        f"Could not redirect to job list after running recipe:\nresponse:{response}")

    def test_recipe_code(self):
        "Test the recipe edit/save code view with POST request"

        pass



class FactoryTest(TestCase):

    def setUp(self):

        owner = models.User.objects.filter(is_superuser=True).first()
        pre = len(models.Project.objects.all())
        self.project = auth.create_project(user=owner, name="test",
                                           text="Text", summary="summary")
        post = len(models.Project.objects.all())
        self.assertTrue(post == (pre + 1), "Error creating project in database")

        self.json_data = {
                "label":"test",
                "help":"Test json data",
                "choices":["test1", "test2"]
        }

        return

    def test_factory_fields(self):
        "Testing factory module that generates fields"

        # All valid field types.
        field_types = factory.get_field_types()

        for display_type in field_types:

            # Test that each field type can be rendered.
            json_data = dict(display=display_type)

            field = factory.dynamic_field(json_data)
            if not field:
                message = f"field generator for display={display_type} failed"
                self.assertFalse(message)

    def test_data_generator(self):
        "Test data generator"

        from biostar.engine import const

        pre = models.Data.objects.count()
        data = auth.create_data(self.project, path=__file__)
        post = models.Data.objects.count()

        self.assertTrue(post == (pre + 1), "Error creating data in database")

        display_type = const.DROPDOWN

        json_data = dict(display=display_type, path=data.get_path())

        field = factory.dynamic_field(json_data, project=self.project)

        if not field:
            message = f"field generator for display={display_type} failed"
            self.assertFalse(message)


class ManagementCommandTest(TestCase):
    def setUp(self):
        self.owner = models.User.objects.filter(is_superuser=True).first()

        pre = len(models.Project.objects.all())
        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary", uid="testing")
        post = len(models.Project.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating project in database")
        self.analysis = auth.create_analysis(project=self.project, json_text='{test:{value:"test"}}',
                                             template="echo {{test.value}}")
        self.job = auth.create_job(analysis=self.analysis)

    def test_add_data(self):
        "Test adding data to a project using management commands "

        pre = models.Data.objects.all().count()
        management.call_command('data', path=__file__, uid="testing")
        post = models.Data.objects.all().count()

        self.assertTrue(post == (pre + 1), "Error creating adding in database with management command")

    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id)

        return

    def Xtest_create_analysis(self):
        "Testing createing an analysis with managment commads"
        pass
