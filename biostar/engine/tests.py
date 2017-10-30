
from django.test import TestCase
from django.test import Client

class SiteTestCase(TestCase):
    def setUp(self):

        pass

    def test_site_pages(self):
        "Checking the site pages."
        c = Client()
        resp = c.post('/')
        self.assertEqual(resp.status_code, 200)



