from django.test import TestCase
from django.test import Client
from django.core.urlresolvers import get_resolver
from .urls import urlpatterns




class SiteTestCase(TestCase):


    def test_site_pages(self):
        "Checking the site pages."
        # Test using demo project and/or admin project
        admin = '1'
        demo = '2'

        idpattern = "(?P<id>\d+)"
        pathpattern = "(?P<path>.+)"

        c = Client()
        #patterns = get_resolver(urlpatterns)
        urls =[]
        for pattern in urlpatterns:
            # check results page.
            #print(pattern)
            pattern = pattern.regex.pattern.replace("^",'').replace("$",'')
            pattern = pattern.replace(idpattern, demo).replace(pathpattern,"results")

            if len(pattern):
                urls.append(pattern)
            else:
                urls.append("/")

        print(urls)
        

        resp = c.post('/')
        self.assertEqual(resp.status_code, 200)



