import logging,os
from django.test import TestCase, override_settings
from unittest.mock import patch, MagicMock
from django.core import management
from django.urls import reverse
from django.conf import settings
from biostar.recipes import auth, const
from biostar.recipes import models, views

from . import util

logger = logging.getLogger('engine')

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'tested'))


@override_settings(MEDIA_ROOT=TEST_ROOT)
class JobViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username=f"tested{util.get_uuid(10)}", email="tested@l.com")
        self.owner.set_password("tested")

        self.project = auth.create_project(user=self.owner, name="tested", text="Text", summary="summary",
                                           uid="tested")

        self.recipe = auth.create_analysis(project=self.project, json_text="{}", template="",
                                           security=models.Analysis.AUTHORIZED)

        self.job = auth.create_job(analysis=self.recipe, user=self.owner)
        self.job.save()



    @patch('biostar.recipes.models.Job.save', MagicMock(name="save"))
    def test_job_edit(self):
        "Test job edit with POST request"

        data = {'name':'tested', 'text':"tested" }
        url  = reverse('job_edit', kwargs=dict(uid=self.job.uid))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.job_edit(request=request, uid=self.job.uid)
        self.process_response(response=response, data=data, save=True)

    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id, verbosity=2)
        management.call_command('job', list=True)

    def test_job_copy_and_paste(self):

        url = reverse("job_copy", kwargs=dict(uid=self.job.uid))
        request = util.fake_request(url=url, data={}, user=self.owner)
        response = views.job_copy(request=request, uid=self.job.uid)
        self.process_response(response, data={})

        board = request.session.get(settings.CLIPBOARD_NAME, {}).get(const.RESULTS_CLIPBOARD, [])
        success = len(board) == 1 and board[0] == self.job.uid
        print(board)

        self.assertTrue(success, "Job uid not copied to clipboard")
        return

    def test_job_delete(self):
        "Test job delete"

        url = reverse('job_delete', kwargs=dict(uid=self.job.uid))

        request = util.fake_request(url=url, data={}, user=self.owner)

        response = views.job_delete(request=request, uid=self.job.uid)

        self.process_response(response=response, data={})

    def test_job_file_copy(self):
        "Test the job file copying interface"

        data = {settings.CLIPBOARD_NAME: const.FILES_CLIPBOARD, 'path': self.job.get_data_dir()}
        copy_url = reverse("job_file_copy", kwargs=dict(uid=self.job.uid, path=self.job.get_data_dir()))
        copy_request = util.fake_request(url=copy_url, data=data, user=self.owner)
        copy_response = views.job_file_copy(request=copy_request, uid=self.job.uid, path=self.job.get_data_dir())
        self.process_response(copy_response, data={})

        paste_url = reverse("file_paste", kwargs=dict(uid=self.project.uid))
        paste_request = util.fake_request(url=paste_url, data=data, user=self.owner)
        paste_response = views.file_paste(request=paste_request, uid=self.project.uid)
        self.process_response(paste_response, data={})

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
                         f"Could not redirect to project view after tested :\nresponse:{response}")

        if save:
            self.assertTrue( models.Job.save.called, "save() method not called")
