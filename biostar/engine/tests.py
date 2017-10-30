from django.test import TestCase
from django.test import Client
from django.core.urlresolvers import get_resolver
from .urls import urlpatterns




class SiteTestCase(TestCase):


    def test_site_pages(self):
        "Checking the site pages."
        admin_project = 1
        demo_project = 2

        idpattern = "(?P<id>\d+)"
        pathpattern = "(?P<path>.+)"


        c = Client()
        #patterns = get_resolver(urlpatterns)
        for pat in urlpatterns:
            print(pat.regex.pattern)
        1/0

        resp = c.post('/')
        self.assertEqual(resp.status_code, 200)



