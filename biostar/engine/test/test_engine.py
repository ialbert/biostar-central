import hjson, logging, os
from django import forms
from django.core import management
from django.test import TestCase
from django.urls import reverse

from biostar.engine import auth, factory
from biostar.engine import models

logger = logging.getLogger('engine')

class ViewsTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = models.User.objects.filter(is_superuser=True).first()

        dbcounter = {
            "project": {"pre": models.Project.objects.count(), "post": models.Project},
            "analysis": {"pre": models.Analysis.objects.count(), "post": models.Analysis},
        }
        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary")
        self.analysis = auth.create_analysis(project=self.project, json_text='{}', template="")

        for model_type, states in dbcounter.items():
            post = states["post"].objects.count()
            self.assertTrue(post == (states["pre"] + 1), f"Error adding {model_type} to database.")

        access = models.Access(user=self.owner,
                      access=models.Access.ADMIN_ACCESS, project=self.project)
        access.save()
        self.project.access_set.add(access)
        self.project.save()


    def test_project_users(self):
        "Test project_users view"

        url = reverse("project_users", kwargs=dict(uid=self.project.uid))
        new_user = models.User.objects.create(email="test@test.com",
                                              first_name="test",
                                              username="test")
        new_user.set_password("test")
        new_user.save()

        info = dict(user=self.owner, project_id=self.project.id,
                    user_id=new_user.id, access=models.Access.ADMIN_ACCESS)

        resp = self.client.post(url, data=info)
        self.assertEqual(resp.status_code, 302)


    def test_project_edit(self):
        "Test project edit"

        url = reverse("project_edit", kwargs=dict(uid=self.project.uid))
        info = dict(text="new text", summary="new summary", name="new name")
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code, 302, f"Error : status_code={resp.status_code} expected 302")


    def test_project_create(self):
        "Test for project creation view"

        info = dict(user=self.owner, name="testing name", summary="test", text="testing")
        resp = self.client.post(reverse("project_create"), info)

        self.assertEqual(resp.status_code, 302, f"Error : status_code={resp.status_code} and expected 302")


    def test_data_edit(self):

        pre = models.Data.objects.count()
        data = auth.create_data(self.project, path=__file__)
        post = models.Data.objects.count()

        self.assertTrue(post == (pre + 1), "Error creating data in database")

        url = reverse("data_edit", kwargs=dict(id=data.id))
        info = dict(summary="new summary", text="new text", name="new name")
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code, 302)


    def test_data_upload(self):
        "Test for data upload interface"

        url = reverse("data_upload", kwargs=dict(uid=self.project.uid))
        info = dict(user=self.owner, summary="test upload", text="test", file=__file__)
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code,302)

    def test_recipe_copy(self):
        new_project = auth.create_project(user=self.owner, name="test",
                                      text="Text", summary="summary")

        url = reverse("recipe_view", kwargs=dict(id=self.analysis.id))
        info = dict(project=new_project.id)
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code, 302)


    def test_recipe_edit(self):
        url = reverse("recipe_edit", kwargs=dict(id=self.analysis.id))
        json_data = {"settings": {"name": "Test"}}
        json_text = hjson.dumps(json_data)

        for option in ("preview", "save"):
            logger.info(f"Testing {option} for recipe_edit")
            info = dict(user=self.owner, text=json_text, save_or_preview=option)
            resp = self.client.post(url, info)

            self.assertEqual(resp.status_code, 302)

    def test_recipe_run(self):
        url = reverse("analysis_run", kwargs=dict(id=self.analysis.id))
        info = dict(user=self.owner)
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code, 302)
        pass

    def Xtest_job_list(self):
        pass

    def Xtest_job_view(self):
        pass

    def Xtest_job_result_view(self):
        pass

    def Xtest_job_file_view(self):
        pass

    def Xtest_job_files_list(self):
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
