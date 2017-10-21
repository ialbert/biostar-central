from django.test import TestCase
from django.test import Client

# 302 is temporary move and 301 is permanent

ALL_PAGES = [("/", 200), ("/projects", 301), ("/login", 301), ("/signup", 301)]


class SimplePageResponses(TestCase):

    def test_pages(self):
        client = Client()

        for url, status in ALL_PAGES:
            response = client.get(url)
            self.assertEqual(response.status_code, status)





