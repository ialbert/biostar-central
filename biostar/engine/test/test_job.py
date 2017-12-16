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

        pre = models.Job.objects.count()
        self.job = auth.create_job(analysis=self.recipe, user=self.owner)
        self.job.save()
        self.assertTrue(models.Job.objects.count() == (pre + 1), "Error creating Job database")


    @patch('biostar.engine.models.Job.save', MagicMock(name="save"))
    def test_job_edit(self):
        "Test job edit with POST request"

        data = {'name':'test', 'text':"testing", 'sticky':True}
        url  = reverse('job_edit', kwargs=dict(id=self.job.id))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.job_edit(request=request, id=self.job.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect after editing job:\nresponse:{response}")

        self.assertTrue(self.job.url() == response.url,
                        f"Could not redirect to correct page: {self.job.url()}!= {response.url}")

        self.assertTrue( models.Job.save.called, "job.save() method not called when editing.")


    def test_job_files_entry(self):
        "Test job_files_entry with POST request"

        management.call_command('job', id=self.job.id)
        url = reverse('job_files_entry', kwargs=dict(id=self.job.id))

        data = {"paths":"run.sh"}

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.job_files_list(request=request, id=self.job.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect after editing job:\nresponse:{response}")
        self.assertTrue(self.job.url() == response.url,
                        f"Could not redirect to correct page: {self.job.url()}!= {response.url}")



    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id, verbosity=2)
        management.call_command('job', list=True)
