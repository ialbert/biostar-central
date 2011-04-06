"""
Main test script executed when we run "manage.py test".
"""
import sys, os, django

if django.VERSION < (1, 3):
    print '*** Django version 1.3 or higher required.'
    print '*** Your version is %s' % '.'.join( map(str, django.VERSION))
    sys.exit()

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

fixture = path(os.path.dirname(__file__), '..', '..', 'home', 'import', 'test-fixture.json' )

from django.test import TestCase
from django.utils import unittest
from django.test.client import Client
from django.conf import settings

class EnvironmentTest(unittest.TestCase):
    # we just need the version to exists
    def test_biostar_version(self):
        self.assertTrue(settings.BIOSTAR_VERSION)

class UrlTest(TestCase):

    fixtures = [ fixture ]

    def test_access(self):
        "Testing that basic URLs function correctly"
        urls = "/ /about/ /member/list/".split()
        c = Client()
        for url in urls:
            resp = c.get(url)
            self.assertEqual(resp.status_code, 200)
    
    def test_redirect(self):
        "Testins redirecting urls"
        urls = "/question/new/".split()
        c = Client()
        for url in urls:
            resp = c.get(url)
            self.assertEqual(resp.status_code, 302)

    def test_login(self):
        "Testing login"
        c = Client()
        r = c.get('/test/login/2/', follow=True)
        self.assertTrue('logout' in r.content)
    
__test__ = { "doctest": """
>>> c = Client()
"""}

