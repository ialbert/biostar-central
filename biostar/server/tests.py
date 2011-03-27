"""
Main test script executed when we run "manage.py test".

>>> c = Client()
>>> resp = c.get('/')
>>> resp.status_code
200
>>> resp = c.get('/about/')
>>> resp.status_code
200
"""
import sys, django
if django.VERSION < (1, 3):
    print '*** Django version 1.3 or higher required.'
    print '*** Your version is %s' % str(django.VERSION)
    sys.exit()

from django.test import TestCase
from django.utils import unittest
from django.test.client import Client

class SimpleTest(TestCase):
    def test_basic_addition(self):
        """
        Tests that 1 + 1 always equals 2.
        """
        self.failUnlessEqual(1 + 1, 2)

__test__ = { "doctest": """
Another way to test that 1 + 1 is equal to 2.

>>> 1 + 1 == 2
True
"""}

