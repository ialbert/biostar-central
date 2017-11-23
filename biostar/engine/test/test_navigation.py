import logging, os
from django.test import TestCase
from django.test import Client
from biostar.engine import auth
from biostar.engine import models

from django.urls import reverse

logger = logging.getLogger('engine')

class SiteNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        user = models.User.objects.all().first()
        self.project = auth.create_project(user=user, name="Test project",
                                           privacy=models.Project.PUBLIC)
        data = auth.create_data(project=self.project, path=__file__)
        analysis = auth.create_analysis(project=self.project, json_text='{}', template="")
        job = auth.create_job(analysis=analysis)

        self.proj_params = dict(id=self.project.id)
        self.analysis_params = dict(id=analysis.id)
        self.data_params = dict(id=data.id)
        self.job_params = dict(id=job.id)

    def visit_urls(self, urls, codes):
        c = Client()
        for url in urls:
            resp = c.get(url)
            if resp.status_code not in codes:
                # print (resp.content)
                # We already know it is an error.
                # Use this to prints the url and the code.
                logger.error(f"")
                logger.error(f"Error accessing: {url}, code={resp.status_code} not in expected values")
                self.assertEqual(url, codes)

    def test_public_pages(self):
        "Checking public pages"

        #TODO:  'job_files_list' not tested yet
        urls = [
            reverse('index'), reverse('info'), reverse('logout'),
            reverse('login'), reverse('signup'),
            reverse('project_list'),
            reverse('data_list', kwargs=self.proj_params),
            reverse('data_view', kwargs=self.data_params),
            reverse('project_view', kwargs=self.proj_params),
            reverse('analysis_list', kwargs= self.proj_params),
            reverse('analysis_view', kwargs=self.analysis_params),
            #reverse('analysis_run', kwargs=self.analysis_params),
            #reverse('analysis_recipe', kwargs=self.analysis_params),
            #reverse('analysis_copy', kwargs=self.analysis_params),
            #reverse('job_list', kwargs=self.proj_params),
            #reverse('job_view', kwargs=self.job_params),
            #reverse('job_files_entry', kwargs=self.job_params),

        ]

        self.visit_urls(urls, [200, 302])

    def test_page_redirect(self):
        "Testing that a redirect occurs for some pages"
        urls = [
            reverse('project_create'),
            reverse('project_edit',  kwargs=self.proj_params),
            reverse('data_upload', kwargs=self.proj_params),
            reverse('data_edit', kwargs=self.data_params),
            reverse('analysis_edit', kwargs=self.analysis_params),
            reverse('job_result_view', kwargs=self.job_params),

        ]

        self.visit_urls(urls, [302])



