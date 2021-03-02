import logging,os
from django.test import TestCase, override_settings
from unittest.mock import patch, MagicMock
from django.core import management
from django.urls import reverse
from django.conf import settings
from biostar.recipes import auth, const
from biostar.recipes import models, views

from biostar.utils.helpers import fake_request, get_uuid

logger = logging.getLogger('engine')

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'tested'))
TOC_ROOT = os.path.join(TEST_ROOT, 'toc')

# Ensure that the table of directory exists.
os.makedirs(TOC_ROOT, exist_ok=True)


@override_settings(MEDIA_ROOT=TEST_ROOT, TOC_ROOT=TOC_ROOT)
class JobViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com")
        self.owner.set_password("tested")

        self.project = auth.create_project(user=self.owner, name="tested", text="Text", summary="summary",
                                           uid="tested")

        self.recipe = auth.create_analysis(project=self.project, json_text="", template="",
                                           security=models.Analysis.AUTHORIZED)

        self.job = auth.create_job(analysis=self.recipe, user=self.owner)
        self.job.save()

    def test_scheduler(self):
        """
        Test task scheduler used to run queued jobs.
        """
        from biostar.recipes.tasks import scheduler
        self.job = auth.create_job(analysis=self.recipe, user=self.owner)
        self.job.state = models.Job.QUEUED
        self.job.save()
        #scheduler

    @patch('biostar.recipes.models.Job.save', MagicMock(name="save"))
    def test_job_edit(self):
        "Test job edit with POST request"

        data = {'name':'tested', 'text':"tested" }
        url  = reverse('job_edit', kwargs=dict(uid=self.job.uid))

        request = fake_request(url=url, data=data, user=self.owner)

        response = views.job_edit(request=request, uid=self.job.uid)
        self.process_response(response=response, data=data, save=True)

    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id, verbosity=2)
        management.call_command('job', list=True)


    def test_job_delete(self):
        "Test job delete"

        url = reverse('job_delete', kwargs=dict(uid=self.job.uid))

        request = fake_request(url=url, data={}, user=self.owner)

        response = views.job_delete(request=request, uid=self.job.uid)

        self.process_response(response=response, data={})

    def test_job_rerun(self):
        "Test Job rerun"
        url = reverse('job_delete', kwargs=dict(uid=self.job.uid))

        request = fake_request(url=url, data={}, user=self.owner)

        response = views.job_rerun(request=request, uid=self.job.uid)
        self.process_response(response=response, data={})


    def test_job_serve(self):
        "Test file serve function."
        from django.http.response import FileResponse

        management.call_command('job', id=self.job.id)

        url = reverse('job_view', kwargs=dict(uid=self.job.uid))

        data = {"paths":"runlog/input.json"}

        request = fake_request(url=url, data=data, user=self.owner)

        response = views.job_serve(request=request, uid=self.job.uid, path=data["paths"])

        self.assertTrue(isinstance(response, FileResponse), "Response is not a file.")


    def process_response(self, response, data, save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to project view after tested :\nresponse:{response}")

        if save:
            self.assertTrue( models.Job.save.called, "save() method not called")
