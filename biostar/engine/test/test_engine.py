import hjson, logging, os
from django import forms
from django.core import management
from django.test import TestCase, RequestFactory
from unittest.mock import patch, MagicMock
from django.urls import reverse
from django.contrib import messages
from biostar.engine import auth, factory
from biostar.engine import models, forms, views

logger = logging.getLogger('engine')



class ProjectTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com")
        self.owner.set_password("test")

        self.factory = RequestFactory()


    @patch('biostar.engine.models.Project.save', MagicMock(name="save"))
    def test_create_view(self):
        "Test project create view"

        # Create fake request
        data = {'name': 'My project', 'uid': 'example'}
        request = self.factory.post(reverse('project_create'), data)
        request.user = self.owner
        request._messages = messages

        # Count before creating one
        pre = models.Project.objects.count()
        response = views.project_create(request)

        print(response, pre, models.Project.objects.count())
        1/0

        pass

    def Xtest_edit_view(self):
        form = YourForm()
        form.cleaned_data = {'name': 'Jan Nowak'}
        form.save()
        pass

    def Xtest_users_view(self):
        form = YourForm()
        form.cleaned_data = {'name': 'Jan Nowak'}
        form.save()
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
