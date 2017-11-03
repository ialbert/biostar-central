
import logging
from django.test import TestCase
from django.test import Client
from biostar.engine import auth
from biostar.engine import models
from biostar.engine import util


from django.urls import reverse

logger = logging.getLogger('engine')


class ProjectTest(TestCase):

    def test_project_create(self):
        "Testing project create"
        owner = models.User.objects.filter(is_superuser=True).first()
        info = dict(user=owner, name="testing name", summary="test",
                    text="testing")

        c = Client()
        resp = c.post(reverse("project_create"), info)

        # Redirects to projects_list
        self.assertEqual(resp.status_code, 302)

    def test_project_edit(self):
        "Testing project editing"

        owner = models.User.objects.filter(is_superuser=True).first()
        test_project = auth.create_project(user=owner, name="test",
                                                text="Text",summary="summary")
        test_project.save()

        new_text = f"new text-{util.get_uuid(4)}"
        new_summary = f"new summary-{util.get_uuid(4)}"
        new_name = f"new name-{util.get_uuid(4)}"
        url = reverse("project_edit", kwargs=dict(id=test_project.id))
        info = dict(text=new_text, summary=new_summary, name=new_name)

        # Emulates user editing here
        c = Client()
        resp = c.post(url, info)

        # Redirects to projects_view
        self.assertEqual(resp.status_code, 302)

    def test_data_upload(self):
        "Test data upload to a project"

        return



class AnalysisTest(TestCase):

    def setUp(self):
        pass

    def test_creation(self):
        pass
    def test_analysis_edit(self):
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