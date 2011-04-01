"""
Main test script executed when we run "manage.py test".
"""
import sys, django

if django.VERSION < (1, 3):
    print '*** Django version 1.3 or higher required.'
    print '*** Your version is %s' % '.'.join( map(str, django.VERSION))
    sys.exit()

from django.test import TestCase
from django.utils import unittest
from django.test.client import Client

class UrlTest(TestCase):
    def test_access(self):
        "Testing that basic URLs function correctly"
        urls = "/ /about/ /members/".split()
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

__test__ = { "doctest": """
Another way to test that 1 + 1 is equal to 2.

>>> 1 + 1 == 2
True
"""}

