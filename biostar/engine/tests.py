import os, logging
from django.test import TestCase
from django.test import Client
from biostar.engine.models import Project, Data, Analysis, Job
from biostar.engine import auth

from django.urls import reverse

logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


class Bunch(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def remove_file(file):
    try:
        os.remove(file)
    except FileNotFoundError:
        pass
    return


def setup_testdata(project):
    file = open("tmp.text", "w")
    file.close()
    data = auth.create_data(project=project, fname="tmp.text")
    remove_file("tmp.text")
    return data


def setup_testanalysis(project):
    analysis = auth.create_analysis(project=project, json_text="{}", template="", user=None, summary='',
                                    name='test', text='')
    analysis.save()
    return analysis


def Xget_urls(current_id, idpattern, pathpattern, data=Bunch(id=1), job=Bunch(id=1),
             analysis=Bunch(id=1), urlpatterns=[]):
    urls = []

    for pattern in urlpatterns:
        pattern = pattern.regex.pattern.replace("^", "/").replace("$", '')
        if "data" in pattern and "data/list" not in pattern:
            current_id = str(data.id)
        elif "analysis" in pattern and "analysis/list" not in pattern:
            current_id = str(analysis.id)
        elif "job" in pattern and "job/list" not in pattern:
            current_id = str(job.id)

        pattern = pattern.replace(idpattern, current_id).replace(pathpattern, "")

        # Ignore admin stuff for now.
        if ("admin" not in pattern):
            urls.append(pattern)

    urls.append("/")
    return urls


class UserNavigation(TestCase):
    def test_foo(self):
        "Testing foo"
        pass

class SiteNavigation(TestCase):

    def setUp(self):
        '''
        logger.setLevel(logging.WARNING)
        self.idpattern = "(?P<id>\d+)"
        self.pathpattern = "(?P<path>.+)/"
        self.project = Project.objects.filter(id=2).first()
        self.data = Data.objects.filter(id=1).first()

        if not self.data:
            self.data = setup_testdata(project=self.project)

        self.analysis = Analysis.objects.filter(id=1).first()

        if not self.analysis:
            self.analysis = setup_testanalysis(project=self.project)

        self.job = Job.objects.filter(id=1).first()

        if not self.job:
            self.job = auth.create_job(analysis=self.analysis)
            self.job.save()
        '''

    def visit_urls(self, urls, code):
        c = Client()
        for url in urls:
            resp = c.get(url)
            if resp.status_code != code:
                # print (resp.content)
                # We already know it is an error.
                # Use this to prints the url and the code.
                self.assertEqual(url, code)

    def test_public_pages(self):
        "Checking public pages"

        params = dict(id=1)

        urls = [
            reverse('index'), reverse('info'),
            reverse('login'), reverse('signup'),
            reverse('project_list'),
            reverse('project_view', kwargs=params),
            reverse('data_list', kwargs=params),
            reverse('analysis_list', kwargs=params),
            reverse('job_list', kwargs=params),
        ]

        self.visit_urls(urls, 200)


    def test_page_redirect(self):
        "Testing that a redirect occurs for some pages"

        params = dict(id=1)
        urls = [
            reverse('data_upload', kwargs=params),
        ]

        self.visit_urls(urls, 300)


    def Xtest_account_pages(self):
        'Tests profile page for user.id=1'

        user_profile_id = '1'

        urls = get_urls(current_id=user_profile_id, idpattern=self.idpattern, pathpattern=self.pathpattern,
                        urlpatterns=accountsurls)

        c = Client()
        accepted = [200, 302]

        accounts_map = {"/accounts/login/": {'email': 'testbuddy@lvh.me', 'password': 'testbuddy@lvh.me'},
                        "/accounts/signup/": {'email': "3@test.com", 'password1': 'test', 'password2': 'test'},
                        "/accounts/logout/": {},
                        f"/accounts/{user_profile_id}/": {},
                        }

        for accounts_page in urls:
            url = f"/accounts{accounts_page}"
            resp = c.post(url, accounts_map.get(url))
            self.assertTrue(resp.status_code in accepted)
            # print(url, resp.status_code)
