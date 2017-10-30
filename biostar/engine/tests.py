from django.test import TestCase
from django.test import Client
from django.core.urlresolvers import get_resolver
from .urls import urlpatterns
from .models import Project, Data, Analysis, Job
from biostar.accounts.urls import urlpatterns as accountsurls



class SiteTestCase(TestCase):

    def setUp(self):
        self.project = Project.objects.filter(id=2).first()
        self.data = Data.project
        pass

    def test_site_pages(self):
        "Checking the site pages."
        # Test using demo project and/or admin project
        demo = '2'

        idpattern = "(?P<id>\d+)"
        pathpattern = "(?P<path>.+)"

        urls =[]
        for pattern in urlpatterns:

            pattern = pattern.regex.pattern.replace("^","/").replace("$",'')

            # Checks results page
            pattern = pattern.replace(idpattern, demo).replace(pathpattern,"results")

            # Ignore admin stuff for now.
            if ("admin" not in pattern):
                urls.append(pattern)

        urls.append("/")
        c = Client()
        for page in urls:

            if page:
                print(page)
                resp = c.post(page)
                print(resp.status_code)
                #self.assertEqual(resp.status_code, 200)



