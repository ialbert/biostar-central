
import logging, os, hjson
from django.test import TestCase
from biostar.engine import auth
from biostar.engine import models
from django.urls import reverse
from django.core import management


logger = logging.getLogger('engine')


class ProjectTest(TestCase):


    def setUp(self):

        logger.setLevel(logging.WARNING)
        self.owner = models.User.objects.filter(is_superuser=True).first()

        pre = len(models.Project.objects.all())

        self.project = auth.create_project(user=self.owner, name="test",
                                                text="Text",summary="summary")
        self.project.save()

        post = len(models.Project.objects.all())

        self.assertTrue(post==(pre+1), "Error creating project in database")

        # using a simple logged in client when needed
        self.client.login(username="1@lvh.me", password="1@lvh.me")


    def test_create_interface(self):
        "Test for project creation interface"

        info = dict(user=self.owner, name="testing name", summary="test",text="testing")
        resp = self.client.post(reverse("project_create"), info, follow=True)

        # Test if user interface works
        self.assertEqual(resp.status_code, 200, f"Error : response.status_code={resp.status_code} and expected 200.")


    def test_edit_interface(self):
        "Test for project editing interface"

        url = reverse("project_edit", kwargs=dict(id=self.project.id))
        info = dict(user=self.owner, text="new text", summary="new summary", name="new name")

        resp = self.client.post(url, info, follow=True)

        self.assertEqual(resp.status_code, 200, f"Error : response.status_code={resp.status_code} and expected 200.")


    def test_data_upload_interface(self):
        "Test for data upload interface"

        url = reverse("data_upload", kwargs=dict(id=self.project.id))
        info = dict(user=self.owner, summary="test upload", text="test", file=__file__)
        resp = self.client.post(url, info, follow=True)

        self.assertEqual(resp.status_code, 200, f"Error : response.status_code={resp.status_code} and expected 200.")


    def test_data_edit_interface(self):
        "Test data edit interface"

        pre = len(models.Data.objects.all())
        data = auth.create_data(self.project, fname=__file__)
        post = len(models.Data.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating data in database")

        url = reverse("data_edit", kwargs=dict(id=data.id))
        info = dict(summary="new summary", text="new text", name="new name")
        resp = self.client.post(url, info, follow=True)

        self.assertEqual(resp.status_code, 200, f"Error : response.status_code={resp.status_code} and expected 200.")


    def test_get_project_list(self):
        "Test getting query of projects belonging to a user ( uses auth.py)"
        projects = auth.get_project_list(user=self.owner)

        return



class DataTest(TestCase):


    def setUp(self):

        self.owner = models.User.objects.filter(is_superuser=True).first()
        pre = len(models.Project.objects.all())
        self.project = auth.create_project(user=self.owner, name="test",
                                                text="Text",summary="summary")
        post = len(models.Project.objects.all())

        self.assertTrue(post==(pre+1), "Error creating project in database")

        pass

    def test_data_create(self):

        "Testing data create with auth"

        pre = len(models.Data.objects.all())
        data = auth.create_data(self.project, fname=__file__)
        post = len(models.Data.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating data in database")
        data_ext= os.path.splitext(data.get_path())[-1]
        file_ext =  os.path.splitext(__file__)[-1]

        self.assertEqual(data_ext, file_ext, f"Created data extension ({data_ext}) does not match the input ({file_ext})")


    def test_data_linkage(self):

        "Test data linkage with auth"

        pre = len(models.Data.objects.all())
        data = auth.create_data(self.project, fname=__file__, link=True)
        post = len(models.Data.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating linked data in database")
        self.assertEqual(data.get_path(), __file__, "Path in database does not math linked file path.")


class AnalysisTest(TestCase):


    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = models.User.objects.filter(is_superuser=True).first()

        dbcounter = {
               "project":{"pre":len(models.Project.objects.all()),"post":models.Project},
               "analysis":{"pre":len(models.Analysis.objects.all()),"post":models.Analysis},
               }
        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary")
        self.analysis = auth.create_analysis(project=self.project, json_text='{}', template="")

        for model_type, states in dbcounter.items():
            post = len(states["post"].objects.all())
            self.assertTrue(post == (states["pre"] + 1), f"Error adding {model_type} to database.")

        # using a simple logged in client when needed
        self.client.login(username="1@lvh.me", password="1@lvh.me")


    def test_analysis_copy_interface(self):
        "Testing analysis copy interface "

        url = reverse("analysis_copy", kwargs=dict(id=self.analysis.id))
        projects = [self.project]
        info = dict(projects=projects)
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code, 200)


    def test_analysis_edit_interface(self):
        "Testing analysis edit interface"

        url = reverse("analysis_edit", kwargs=dict(id=self.analysis.id))
        json_data = {"settings":{"name":"Test"}}
        json_text = hjson.dumps(json_data)

        for option in ("preview", "save"):

            logger.info(f"Testing {option} for analysis_edit")
            info = dict(user=self.owner, text=json_text, save_or_preview=option)
            resp = self.client.post(url, info, follow=True)

            self.assertEqual(resp.status_code, 200)

    def test_analysis_run_interface(self):
        "Testing analysis run interface"

        url = reverse("analysis_run", kwargs=dict(id=self.analysis.id))
        info=dict(user=self.owner)
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code, 200)


class JobTest(TestCase):

    def setUp(self):

        logger.setLevel(logging.WARNING)

        self.owner = models.User.objects.filter(is_superuser=True).first()

        dbcounter = {
               "project":{"pre":len(models.Project.objects.all()),"post":models.Project},
               "analysis":{"pre":len(models.Analysis.objects.all()),"post":models.Analysis},
               "job":{"pre":len(models.Job.objects.all()),"post":models.Job}
               }

        self.project = auth.create_project(user=self.owner, name="test",text="Text", summary="summary")
        self.analysis = auth.create_analysis(project=self.project, json_text='{test:{value:"test"}}',
                                             template="echo {{test.value}}")
        self.job = auth.create_job(analysis=self.analysis)

        for model_type, states in dbcounter.items():
            post = len(states["post"].objects.all())
            self.assertTrue( post == (states["pre"]+1), f"Error adding {model_type} to database." )

        # using a simple logged in client when needed
        self.client.login(username="1@lvh.me", password="1@lvh.me")

    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id)

        return



class TasksTests(TestCase):

    def setUp(self):
        return



class CommandTests(TestCase):


    def setUp(self):

        self.owner = models.User.objects.filter(is_superuser=True).first()

        pre = len(models.Project.objects.all())
        self.project = auth.create_project(user=self.owner, name="test",
                                                text="Text",summary="summary", uid="testing")
        post = len(models.Project.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating project in database")



    def test_add_data(self):

        "Test adding data to a project using management commands "
        pre = len(models.Data.objects.all())
        management.call_command('data', path = __file__, uid="testing")
        post = len(models.Data.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating adding in database with management command")
        pass

    def test_link_data(self):
        "Test linking data to a project using management commands "
        pre = len(models.Data.objects.all())

        management.call_command('data', path = __file__, uid="testing", link=True)

        post = len(models.Data.objects.all())

        self.assertTrue(post == (pre + 1), "Error creating adding in database with management command")


    def test_analysis_add(self):
        "Test adding analysis to project using management commands"

        pass

    def test_create_project(self):
        "Test creating project using a json management command"

        # make a temp
        #fp = tempfile.TemporaryFile()
        #write a json file into tmp file and get it done.
        #fp.write(stream.read(CHUNK))
        #fp.seek(0)
        #stream = File(fp)

        pass






