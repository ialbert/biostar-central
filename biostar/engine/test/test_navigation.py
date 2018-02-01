import logging, os
from django.test import TestCase
from django.test import Client
from biostar.engine import auth
from biostar.engine import models

from django.urls import reverse
from . import util

logger = logging.getLogger('engine')

class SiteNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        self.owner = models.User.objects.create(username="test", email="test@test.com")
        self.owner.set_password("testing")
        self.owner.save()

        self.project = auth.create_project(user=self.owner, name="Test project",
                                           privacy=models.Project.PUBLIC, uid="testing")
        data = auth.create_data(project=self.project, path=__file__)
        analysis = auth.create_analysis(project=self.project, json_text='{}', template="")
        self.job = auth.create_job(analysis=analysis)

        self.proj_params = dict(uid=self.project.uid)
        self.analysis_params = dict(uid=analysis.uid)
        self.data_params = dict(uid=data.uid)
        self.job_params = dict(uid=self.job.uid)

    def visit_urls(self, urls, codes):
        c = Client()
        c.login(username="test", email='test@test.com', password='testing')
        for url in urls:

            resp = c.get(url, data={"q":"test"})
            if resp.status_code not in codes:
                # print (resp.content)
                # We already know it is an error.
                # Use this to prints the url and the code.
                logger.error(f"")
                logger.error(f"Error accessing: {url}, code={resp.status_code} not in expected values")
                self.assertEqual(url, codes)

    def test_public_pages(self):
        "Checking public pages"

        urls = [
            reverse('index'),
            reverse('docs', kwargs=dict(name='info')),
            reverse('logout'),
            reverse('login'),
            reverse('project_list'),
            reverse('data_list', kwargs=self.proj_params),
            reverse('data_view', kwargs=self.data_params),
            reverse("data_nav", kwargs=self.proj_params),
            reverse("data_files_entry", kwargs=self.data_params),
            reverse('project_view', kwargs=self.proj_params),
            reverse('project_users', kwargs=self.proj_params),
            reverse('project_create'),
            reverse('project_edit', kwargs=self.proj_params),
            reverse('recipe_list', kwargs= self.proj_params),
            reverse('recipe_view', kwargs=self.analysis_params),
            reverse("recipe_code", kwargs=self.analysis_params),
            reverse('recipe_create', kwargs=self.proj_params),
            reverse('recipe_run', kwargs=self.analysis_params),
            reverse('recipe_view', kwargs=self.analysis_params),
            reverse('job_list', kwargs=self.proj_params),
            reverse('job_view', kwargs=self.job_params),
            reverse('job_edit', kwargs=self.job_params),
            reverse('job_files_entry', kwargs=self.job_params),
            reverse('data_upload', kwargs=self.proj_params),
            reverse('data_edit', kwargs=self.data_params),
            reverse('recipe_edit', kwargs=self.analysis_params),

        ]

        self.visit_urls(urls, [200])

    def test_page_redirect(self):
        "Testing that a redirect occurs for some pages"
        urls = [
            reverse('signup'),
            reverse("data_copy", kwargs=self.data_params),
            reverse("recipe_copy", kwargs=self.analysis_params),
            reverse("data_paste", kwargs=self.proj_params),
        ]

        self.visit_urls(urls, [302, 200])



