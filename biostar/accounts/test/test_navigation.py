import logging, os
from django.test import Client
from django.test import TestCase
from biostar.engine import models
from biostar.engine import auth

logger = logging.getLogger('engine')


class AccountsNavigation(TestCase):

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

        data = auth.create_data(project=self.project, path=__file__)
        analysis = auth.create_analysis(project=self.project, json_text='{}', template="")
        job = auth.create_job(analysis=analysis)

        proj_params = dict(id=self.project.id)
        analysis_params = dict(id=analysis.id)
        data_params = dict(id=data.id)
        job_params = dict(id=job.id)
