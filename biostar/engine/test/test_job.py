import logging,os
from django.test import TestCase, override_settings
from unittest.mock import patch, MagicMock
from django.core import management
from django.urls import reverse
from django.conf import settings
from biostar.engine import auth
from biostar.engine import models, views

from . import util

logger = logging.getLogger('engine')

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, '..', 'export', 'test'))


@override_settings(MEDIA_ROOT=TEST_ROOT)
class JobViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com")
        self.owner.set_password("test")

        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary",
                                           uid="testing")

        self.recipe = auth.create_analysis(project=self.project, json_text="{}", template="",
                                           security=models.Analysis.AUTHORIZED)

        self.job = auth.create_job(analysis=self.recipe, user=self.owner)
        self.job.save()


    @patch('biostar.engine.models.Job.save', MagicMock(name="save"))
    def test_job_edit(self):
        "Test job edit with POST request"

        data = {'name':'test', 'text':"testing", 'sticky':True}
        url  = reverse('job_edit', kwargs=dict(uid=self.job.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.job_edit(request=request, uid=self.job.uid)
        self.process_response(response=response, data=data, save=True)


    def test_job_files_copy(self):
        "Test files copy with POST request"

        management.call_command('job', id=self.job.id)
        url = reverse('job_files_entry', kwargs=dict(uid=self.job.uid))

        data = {"paths":"runlog/input.json"}

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.job_files_list(request=request, uid=self.job.uid)

        self.process_response(response=response, data=data)

    @patch('biostar.engine.models.Data.save', MagicMock(name="save"))
    def test_job_files_paste(self):

        management.call_command('job', id=self.job.id)

        url = reverse("files_paste", kwargs=dict(uid=self.project.uid))

        data = {}
        request = util.fake_request(url=url, data=data, user=self.owner)
        request.session["files_clipboard"] = [auth.join(self.job.path, "runlog")]

        request.session["files_clipboard"].append(self.job.uid)

        response = views.files_paste(request=request, uid=self.project.uid)

        self.process_response(response=response, data=data)

        self.assertTrue(models.Data.save.called, "save() method not called when editing.")


    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id, verbosity=2)
        management.call_command('job', list=True)


    def process_response(self, response, data, save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to project view after testing :\nresponse:{response}")

        if save:
            self.assertTrue( models.Job.save.called, "save() method not called")