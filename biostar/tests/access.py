"""
Main test script is executed when running::

    biostar.sh test"

"""
import sys, os, django

from django.test import TestCase
from django.utils import unittest
from django.test.client import Client
from django.conf import settings

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

fixture = path(os.path.dirname(__file__), '..', '..', 'home', 'import', 'test-fixture.json' )

class InternalTest(TestCase):
    "Tests the views access from within Django"

    fixtures = [ fixture ]

    def test_biostar_version(self):
        self.assertTrue(settings.BIOSTAR_VERSION)

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

def suite():
    s = unittest.TestLoader().loadTestsFromTestCase(InternalTest)
    return s
