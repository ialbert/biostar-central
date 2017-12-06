import hjson, logging, os
from django import forms
from django.core import management
from django.test import TestCase
from django.urls import reverse

from biostar.engine import auth, factory
from biostar.engine import models

logger = logging.getLogger('engine')


class ProjectTest(TestCase):
    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = models.User.objects.filter(is_superuser=True).first()

        pre = len(models.Project.objects.all())

        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary")
        self.project.save()

        post = len(models.Project.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating project in database")

        # using a simple logged in client when needed
        self.client.login(username="1@lvh.me", password="1@lvh.me")


class DataTest(TestCase):
    def setUp(self):
        self.owner = models.User.objects.filter(is_superuser=True).first()
        pre = len(models.Project.objects.all())
        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary")
        post = len(models.Project.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating project in database")

        pass

    def test_data_linkage(self):
        "Test data linkage with auth"

        pre = len(models.Data.objects.all())
        data = auth.create_data(self.project, path=__file__)
        post = len(models.Data.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating linked data in database")


class AnalysisTest(TestCase):
    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = models.User.objects.filter(is_superuser=True).first()

        dbcounter = {
            "project": {"pre": len(models.Project.objects.all()), "post": models.Project},
            "analysis": {"pre": len(models.Analysis.objects.all()), "post": models.Analysis},
        }
        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary")
        self.analysis = auth.create_analysis(project=self.project, json_text='{}', template="")

        for model_type, states in dbcounter.items():
            post = len(states["post"].objects.all())
            self.assertTrue(post == (states["pre"] + 1), f"Error adding {model_type} to database.")

        # using a simple logged in client when needed
        self.client.login(username="1@lvh.me", password="1@lvh.me")

    def test_recipe_copy_interface(self):
        "Testing analysis copy interface "
        new_project = auth.create_project(user=self.owner, name="test",
                                      text="Text", summary="summary")

        url = reverse("recipe_copy", kwargs=dict(id=self.analysis.id))
        info = dict(project=new_project.id)
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code, 302)

        #1/0

    def test_recipe_edit_interface(self):
        "Testing analysis edit interface"

        url = reverse("recipe_edit", kwargs=dict(id=self.analysis.id))
        json_data = {"settings": {"name": "Test"}}
        json_text = hjson.dumps(json_data)

        for option in ("preview", "save"):
            logger.info(f"Testing {option} for recipe_edit")
            info = dict(user=self.owner, text=json_text, save_or_preview=option)
            resp = self.client.post(url, info)

            self.assertEqual(resp.status_code, 302)

    def Xtest_analysis_run_interface(self):
        "Testing analysis run interface"

        url = reverse("analysis_run", kwargs=dict(id=self.analysis.id))
        info = dict(user=self.owner)
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code, 302)


class JobTest(TestCase):
    def setUp(self):
        logger.setLevel(logging.WARNING)

        self.owner = models.User.objects.filter(is_superuser=True).first()

        dbcounter = {
            "project": {"pre": len(models.Project.objects.all()), "post": models.Project},
            "analysis": {"pre": len(models.Analysis.objects.all()), "post": models.Analysis},
            "job": {"pre": len(models.Job.objects.all()), "post": models.Job}
        }

        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary")
        self.analysis = auth.create_analysis(project=self.project, json_text='{test:{value:"test"}}',
                                             template="echo {{test.value}}")
        self.job = auth.create_job(analysis=self.analysis)

        for model_type, states in dbcounter.items():
            post = len(states["post"].objects.all())
            self.assertTrue(post == (states["pre"] + 1), f"Error adding {model_type} to database.")

    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id)

        return


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

        pre = len(models.Data.objects.all())
        data = auth.create_data(self.project, path=__file__)
        post = len(models.Data.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating data in database")

        display_type = const.DROPDOWN

        json_data = dict(display=display_type, path=data.get_path())

        field = factory.dynamic_field(json_data, project=self.project)

        if not field:
            message = f"field generator for display={display_type} failed"
            self.assertFalse(message)



class CommandTests(TestCase):
    def setUp(self):
        self.owner = models.User.objects.filter(is_superuser=True).first()

        pre = len(models.Project.objects.all())
        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary", uid="testing")
        post = len(models.Project.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating project in database")

    def test_add_data(self):
        "Test adding data to a project using management commands "

        pre = models.Data.objects.all().count()
        management.call_command('data', path=__file__, uid="testing")
        post = models.Data.objects.all().count()

        self.assertTrue(post == (pre + 1), "Error creating adding in database with management command")


    def test_link_data(self):
        "Test linking data to a project using management commands "
        pre = len(models.Data.objects.all())

        management.call_command('data', path=__file__, uid="testing")

        post = len(models.Data.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating adding in database with management command")


class ViewsTest(TestCase):

    def setUp(self):

        self.owner = models.User.objects.filter(is_superuser=True).first()
        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary")
        self.project.save()


    def test_project_users(self):
        "Test project_users view"

        url = reverse("project_users", kwargs=dict(uid=self.project.uid))
        new_user = models.User.objects.create(email="test@test.com",
                                              first_name="test",
                                              username="test")
        new_user.set_password("test")
        new_user.save()

        info = dict(user=self.owner, add_or_remove="add", users=new_user.id)

        resp = self.client.post(url, data=info)

        #1/0
        pass


    def test_project_edit(self):
        "Test project edit"
        url = reverse("project_edit", kwargs=dict(uid=self.project.uid))
        info = dict(text="new text", summary="new summary", name="new name")

        resp = self.client.post(url, info, follow=True)

        #self.assertEqual(resp.status_code, 200, f"Error : response.status_code={resp.status_code} and expected 200.")
        pass

    def test_project_create(self):
        "Test for project creation view"

        info = dict(user=self.owner, name="testing name", summary="test", text="testing")
        resp = self.client.post(reverse("project_create"), info, follow=True)

        # Test if user interface works
        #self.assertEqual(resp.status_code, 200, f"Error : response.status_code={resp.status_code} and expected 200.")
        pass

    def test_data_edit(self):

        pre = len(models.Data.objects.all())
        data = auth.create_data(self.project, path=__file__)
        post = len(models.Data.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating data in database")

        url = reverse("data_edit", kwargs=dict(id=data.id))
        info = dict(summary="new summary", text="new text", name="new name")
        resp = self.client.post(url, info, follow=True)

        #self.assertEqual(resp.status_code, 200, f"Error : response.status_code={resp.status_code} and expected 200.")

        pass

    def test_data_upload(self):
        "Test for data upload interface"

        url = reverse("data_upload", kwargs=dict(uid=self.project.uid))
        info = dict(user=self.owner, summary="test upload", text="test", file=__file__)
        resp = self.client.post(url, info, follow=True)

       # self.assertEqual(resp.status_code, 200, f"Error : response.status_code={resp.status_code} and expected 200.")

    def Xtest_analysis_list(self):
        pass


    def Xtest_recipe_view(self):
        pass


    def Xtest_analysis_recipe(self):
        pass


    def Xtest_analysis_copy(self):
        pass

    def Xtest_analysis_run(self):
        pass


    def Xtest_recipe_edit(self):
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







