import logging, re

from django.test import TestCase, Client
from django.core.urlresolvers import reverse
from django.conf import settings

VERBOSE = True

class ClientTestCase(TestCase):
    HOST = "www.lvh.me:8080"

    haystack = logging.getLogger("haystack")

    def setUp(self):
        # Need to disable the haystack logger
        self.haystack.setLevel(logging.CRITICAL)
        self._client = Client(HTTP_HOST=ClientTestCase.HOST)

    def tearDown(self):
        self.haystack.setLevel(logging.WARNING)

    def get(self, url, code=200, kwargs={}, data={}, follow=False, pattern=None):
        r = self._client.get(reverse(url, kwargs=kwargs), data, follow=follow)
        self.assertEqual(r.status_code, code)
        if pattern:
            content = r.content.decode("utf-8")
            found  = re.search(pattern, content, re.IGNORECASE)
            if not found and VERBOSE:
                print(pattern)
                print ("---")
                print (content)
                print ("---")
            self.assertTrue(found)
        return r

    def post(self, url, kwargs={}, data={}, pattern=None, follow=False):
        r = self._client.post(reverse(url, kwargs=kwargs), data, follow=follow)

        if follow:
            self.assertEqual(r.status_code, 200)
        else:
            self.assertEqual(r.status_code, 302)

        if pattern:
            content = r.content.decode("utf-8")
            result = re.search(pattern, content, re.IGNORECASE)
            if not result:
                print (content)
                print("*** unable to find pattern: {0} in content.".format(pattern))
                self.assertTrue(result)

        return r