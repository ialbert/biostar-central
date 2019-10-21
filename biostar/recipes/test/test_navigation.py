import logging, os
from django.test import TestCase
from django.test import Client
from biostar.recipes import auth
from biostar.recipes import models
from django.conf import settings

from . import util
from django.urls import reverse


logger = logging.getLogger('engine')


class SiteNavigation(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.username = f"tested{util.get_uuid(10)}"

        self.owner = models.User.objects.create(username=self.username, is_staff=True,
                                                email="tested@tested.com")
        self.owner.set_password("tested")
        self.owner.save()

        self.project = auth.create_project(user=self.owner, name="Test project",
                                           privacy=models.Project.PUBLIC, uid="tested")
        data = auth.create_data(project=self.project, path=__file__)
        analysis = auth.create_analysis(project=self.project, json_text='{}', template="")
        self.job = auth.create_job(analysis=analysis)

        self.proj_params = dict(uid=self.project.uid)
        self.analysis_params = dict(uid=analysis.uid)
        self.data_params = dict(uid=data.uid)
        self.job_params = dict(uid=self.job.uid)

    def visit_urls(self, urls, codes, anon_urls=[]):
        c = Client()
        # Used to test norname urls
        c.login(username=self.username, email='tested@tested.com', password='tested')
        # Used to test with anon users
        anon_c = Client()

        def visit(pages, client):
            for url in pages:
                print(url)
                resp = client.get(url, data={"q": "tested"})
                code = resp.status_code
                if code not in codes:
                    # We already know it is an error.
                    # Use this to prints the url and the code.
                    logger.error(f"")
                    logger.error(f"Error accessing: {url}, code={code} not in expected values {codes}")
                    self.assertTrue(code in codes)

        visit(pages=urls, client=c)
        visit(pages=anon_urls, client=anon_c)

    def test_public_pages(self):
        "Checking public pages"

        ajax_urls = [
            reverse('add_vars'),
            reverse('recipe_fields', kwargs=self.analysis_params),
            reverse('preview_template'),
            reverse('preview_json')
        ]

        api_urls = [

            reverse('api_list'),
            reverse('recipe_api_json', kwargs=self.analysis_params),
            reverse('recipe_api_template', kwargs=self.analysis_params)
        ]
        anon_urls = [
            reverse("index"),
            reverse("project_list"),
            reverse('project_view', kwargs=self.proj_params)
        ]

        urls = [
            reverse('index'),
            reverse('logout'),
            reverse('login'),
            reverse('search'),
            reverse('project_list'),
            reverse('project_list_private'),
            reverse('data_list', kwargs=self.proj_params),
            reverse('data_view', kwargs=self.data_params),
            reverse('data_upload', kwargs=self.proj_params),
            reverse('data_edit', kwargs=self.data_params),
            reverse('project_view', kwargs=self.proj_params),
            reverse('project_users', kwargs=self.proj_params),
            reverse('project_info', kwargs=self.proj_params),
            reverse('project_create'),
            reverse('project_edit', kwargs=self.proj_params),
            reverse('recipe_list', kwargs=self.proj_params),
            reverse('recipe_view', kwargs=self.analysis_params),
            reverse('recipe_run', kwargs=self.analysis_params),
            reverse('recipe_view', kwargs=self.analysis_params),
            reverse('recipe_edit', kwargs=self.analysis_params),
            reverse('job_list', kwargs=self.proj_params),
            reverse('job_view', kwargs=self.job_params),
            reverse('job_edit', kwargs=self.job_params),

        ]

        self.visit_urls(urls=urls, codes=[200])
        self.visit_urls(urls=api_urls, codes=[200])
        self.visit_urls(urls=ajax_urls, codes=[200])
        self.visit_urls(anon_urls=anon_urls, urls=[], codes=[200])

    def test_page_redirect(self):
        "Testing that a redirect occurs for some pages"
        urls = [
            reverse('signup'),
            reverse("recycle_bin"),
        ]

        self.visit_urls(urls, [302, 200])



