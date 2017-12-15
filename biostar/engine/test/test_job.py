import logging
from django.test import TestCase
from unittest.mock import patch, MagicMock
from django.core import management
from django.urls import reverse

from biostar.engine import auth
from biostar.engine import models, views

from . import util

logger = logging.getLogger('engine')


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
        self.assertTrue(models.Job.objects.count() == (pre + 1), "Error creating Job database")


    @patch('biostar.engine.models.Job.save', MagicMock(name="save"))
    def test_job_edit(self):
        "Test job edit with POST request"

        data = {'name':'test', 'text':"testing", 'sticky':True}
        url  = reverse('job_edit', kwargs=dict(id=self.job.id))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.job_edit(request=request, id=self.job.id)

        util.remove_test_folders(self.project.get_project_dir())
        util.remove_test_folders(self.job.path)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect after editing job:\nresponse:{response}")

        self.assertTrue(self.job.url() == response.url,
                        f"Could not redirect to correct page: {self.job.url()}!= {response.url}")

        self.assertTrue( models.Job.save.called, "job.save() method not called when editing.")


    def test_job_files_entry(self):
        "Test job_files_entry with POST request"

        util.remove_test_folders(self.project.get_project_dir())
        util.remove_test_folders(self.job.path)

        # Create files in job dir to copy
        #management.call_command('job', id=self.job.id)

        pass

    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id, verbosity=2)
        management.call_command('job', list=True)

        util.remove_test_folders(self.project.get_project_dir())
        util.remove_test_folders(self.job.path)