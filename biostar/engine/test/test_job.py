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
    def test_job_toggle(self):
        "Test object_state_toggle using a job. "

        url = reverse('toggle_state', kwargs=dict(uid=self.job.uid, obj_type='job'))

        request = util.fake_request(url=url, data={}, user=self.owner)

        response = views.object_state_toggle(request=request, uid=self.job.uid, obj_type="job")
        self.process_response(response=response, data={}, save=True)

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
        url = reverse('job_view', kwargs=dict(uid=self.job.uid))

        data = {"paths":"runlog/input.json"}

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.job_view(request=request, uid=self.job.uid)

        self.process_response(response=response, data=data)

        # Test clear clipboard view
        views.clear_clipboard(request=request, uid=self.project.uid, board="files_clipboard")

        self.assertTrue(request.session.get("files_clipboard")==None, "Clear clipboard not working")


    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id, verbosity=2)
        management.call_command('job', list=True)


    def test_job_serve(self):
        "Test file serve function."
        from django.http.response import FileResponse

        management.call_command('job', id=self.job.id)

        url = reverse('job_view', kwargs=dict(uid=self.job.uid))

        data = {"paths":"runlog/input.json"}

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.job_serve(request=request, uid=self.job.uid, path=data["paths"])

        self.assertTrue(isinstance(response, FileResponse), "Response is not a file.")


    def process_response(self, response, data, save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to project view after testing :\nresponse:{response}")

        if save:
            self.assertTrue( models.Job.save.called, "save() method not called")