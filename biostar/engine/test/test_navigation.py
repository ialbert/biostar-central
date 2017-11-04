import logging, os
from django.test import TestCase
from django.test import Client
from biostar.engine import auth
from biostar.engine import models

from django.urls import reverse

logger = logging.getLogger('engine')


class UserAccountTests(TestCase):

    def test_pages(self):
        'Testing login,logout, and signup'
        urls = [
            (reverse('signup'), dict(info=dict(email="test1@lvh.me", password="test1@lvh.me"),code=200)),
            (reverse('login'),dict(info=dict(email="test1@lvh.me", password="test1@lvh.me"),code=200)),
            (reverse('logout'), dict(info=dict(), code=302)),
        ]

        c = Client()

        for url, info in urls:
            resp = c.post(url, info["info"])
            self.assertEqual(resp.status_code, info["code"])



class SiteNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        user = models.User.objects.all().first()
        self.project = auth.create_project(user=user, name="Test project")

    def visit_urls(self, urls, code):
        c = Client()
        for url in urls:
            resp = c.get(url)
            if resp.status_code != code:
                # print (resp.content)
                # We already know it is an error.
                # Use this to prints the url and the code.
                logger.error(f"")
                logger.error(f"Error accessing: {url}, code={resp.status_code}")
                self.assertEqual(url, code)

    def test_public_pages(self):
        "Checking public pages"

        data = auth.create_data(project=self.project, fname=__file__)
        analysis = auth.create_analysis(project=self.project, json_text='{}', template="")
        job = auth.create_job(analysis=analysis)

        proj_params = dict(id=self.project.id)
        analysis_params = dict(id=analysis.id)
        data_params = dict(id=data.id)
        job_params = dict(id=job.id)

        #TODO:  'job_files_list' not tested yet
        urls = [
            reverse('index'), reverse('info'), reverse('logout'),
            reverse('login'), reverse('signup'),
            reverse('project_list'),
            reverse('project_view', kwargs=proj_params),
            reverse('data_list', kwargs=proj_params),
            reverse('data_view', kwargs=data_params),
            reverse('analysis_list', kwargs=proj_params),
            reverse('analysis_view', kwargs=analysis_params),
            reverse('analysis_run', kwargs=analysis_params),
            reverse('analysis_recipe', kwargs=analysis_params),
            reverse('analysis_copy', kwargs=analysis_params),
            reverse('job_list', kwargs=proj_params),
            reverse('job_view', kwargs=job_params),
            reverse('job_files_entry', kwargs=job_params),
            reverse('job_result_view', kwargs=job_params),
        ]

        self.visit_urls(urls, 200)


    def test_page_redirect(self):
        "Testing that a redirect occurs for some pages"
        params = dict(id=1)
        urls = [
            reverse('project_create'),
            reverse('project_edit',  kwargs=params),
            reverse('data_upload', kwargs=params),
            reverse('data_edit', kwargs=params),
            reverse('analysis_edit', kwargs=params),
        ]

        self.visit_urls(urls, 302)



