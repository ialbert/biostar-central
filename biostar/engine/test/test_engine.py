
import logging, os, hjson
from django.test import TestCase
from biostar.engine import auth
from biostar.engine import models
from django.urls import reverse


logger = logging.getLogger('engine')


class ProjectTest(TestCase):

    def setUp(self):

        self.owner = models.User.objects.filter(is_superuser=True).first()
        self.project = auth.create_project(user=self.owner, name="test",
                                                text="Text",summary="summary")
        self.project.save()

        # using a simple logged in client when needed
        self.client.login(username="1@lvh.me", password="1@lvh.me")


    def test_project_create(self):
        "Testing project create"

        # Create with no image
        info = dict(user=self.owner, name="testing name", summary="test",text="testing")
        resp = self.client.post(reverse("project_create"), info)

        # Redirects to projects_list
        self.assertEqual(resp.status_code, 302)


    def test_project_edit(self):
        "Testing project editing"

        url = reverse("project_edit", kwargs=dict(id=self.project.id))
        info = dict(user=self.owner, text="new text", summary="new summary", name="new name")

        resp = self.client.post(url, info, follow=True)

        self.assertEqual(resp.status_code, 200)


    def test_data_upload(self):
        "Test data upload form to a sample project"

        url = reverse("data_upload", kwargs=dict(id=self.project.id))
        info = dict(user=self.owner, summary="test upload", text="test", file=__file__)
        resp = self.client.post(url, info, follow=True)

        self.assertEqual(resp.status_code, 200)


    def test_data_edit(self):
        "Test data edit"

        data = auth.create_data(self.project, fname=__file__)

        url = reverse("data_edit", kwargs=dict(id=data.id))
        info = dict(summary="new summary", text="new text", name="new name")
        resp = self.client.post(url, info, follow=True)

        self.assertEqual(resp.status_code, 200)



class DataTest(TestCase):

    def setUp(self):
        pass

    def test_unpack(self):
        "Testing data unpack using tasks"
        return

    def test_copy(self):
        "Testing data copy using tasks"
        return



class AnalysisTest(TestCase):


    def setUp(self):
        self.owner = models.User.objects.filter(is_superuser=True).first()
        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary")
        self.project.save()
        self.analysis = auth.create_analysis(project=self.project, json_text='{}', template="")
        self.analysis.save()

        # using a simple logged in client when needed
        self.client.login(username="1@lvh.me", password="1@lvh.me")

    #def test_analysis_create(self):
    #    "Testing analysis creation form ( not implemented yet)  "
    #    pass


    def test_analysis_copy(self):
        "Testing analysis copy"

        url = reverse("analysis_copy", kwargs=dict(id=self.analysis.id))
        projects = [self.project]
        info = dict(projects=projects)
        resp = self.client.post(url, info)

        self.assertEqual(resp.status_code, 200)


    def test_analysis_edit(self):
        "Testing analysis edit"

        url = reverse("analysis_edit", kwargs=dict(id=self.analysis.id))
        json_data = {"settings":{"name":"Test"}}
        json_text = hjson.dumps(json_data)

        for option in ("preview", "save"):

            logger.info(f"Testing {option} for analysis_edit")
            info = dict(user=self.owner, text=json_text, save_or_preview=option)
            resp = self.client.post(url, info, follow=True)

            self.assertEqual(resp.status_code, 200)


    def test_analysis_run(self):
        "Testing analysis run"

        pass



class JobTest(TestCase):

    def setUp(self):
        pass

    def test_creation(self):
        "Test Job creation"
        return

    def test_results(self):
        "Test results page"
        return