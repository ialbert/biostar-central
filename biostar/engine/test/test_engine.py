
import logging, os
from django.test import TestCase
from biostar.engine import auth
from biostar.engine import models
from django.urls import reverse


logger = logging.getLogger('engine')


class ProjectTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = models.User.objects.filter(is_superuser=True).first()
        self.project = auth.create_project(user=self.owner, name="test",
                                                text="Text",summary="summary")
        self.project.save()

        # using a simple logged in client when needed
        self.client.login(username="1@lvh.me", password="1@lvh.me")


    def test_project_create(self):
        "Testing project create form"

        # Create with no image
        info = dict(user=self.owner, name="testing name", summary="test",text="testing")
        resp = self.client.post(reverse("project_create"), info)

        # Redirects to projects_list
        self.assertEqual(resp.status_code, 302)


    def test_project_edit(self):
        "Testing project editing form"

        url = reverse("project_edit", kwargs=dict(id=self.project.id))
        info = dict(user=self.owner, text="new text", summary="new summary", name="new name")

        resp = self.client.post(url, info, follow=True)

        self.assertEqual(resp.status_code, 200)


    def test_data_upload(self):
        "Test data upload form to a sample project"

        url = reverse("data_upload", kwargs=dict(id=self.project.id))
        test_file= open("test", "w")
        test_file.close()
        os.remove("test")

        info = dict(user=self.owner, summary="test upload", text="test", file="test")
        resp = self.client.post(url, info, follow=True)

        self.assertEqual(resp.status_code, 200)


    def test_data_edit(self):
        "Test data edit form "

        data = auth.create_data(self.project, fname=__file__)

        url = reverse("data_edit", kwargs=dict(id=data.id))
        info = dict(summary="new summary", text="new text", name="new name")
        resp = self.client.post(url, info, follow=True)

        self.assertEqual(resp.status_code, 200)


class DataTest(TestCase):

    def setUp(self):
        pass

    def test_unpack(self):
        "Test data unpack using tasks"
        pass

    def test_copy(self):
        "Test data copy using tasks"



class AnalysisTest(TestCase):

    def setUp(self):
        pass

    def test_creation(self):
        "Test analysis creation "
        pass

    def test_analysis_edit(self):
        "Test analysis edit"
        pass

    def test_analysis_run(self):
        "Test analysis run"
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