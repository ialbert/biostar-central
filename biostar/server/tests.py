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


# Twill based functional tests
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import md5, sha

import twill
from twill import commands as tc
from django.core.servers.basehttp import AdminMediaHandler
from django.core.handlers.wsgi import WSGIHandler
from StringIO import StringIO

def twill_setup():
    app = AdminMediaHandler(WSGIHandler())
    twill.add_wsgi_intercept(DOMAIN, PORT, lambda: app)

def twill_teardown():
    twill.remove_wsgi_intercept(DOMAIN, PORT)

# where to run the test server
DOMAIN, PORT = '127.0.0.1', 8080
HOME_PAGE = "http://%s:%s" % (DOMAIN, PORT)

class FunctionalTest(TestCase):
    
    fixtures = [ fixture ]

    def setUp(self):
        twill_setup()
        
    def tearDown(self):
        twill_teardown()
    
    def test_home(self):

        # go to the home page
        tc.go(HOME_PAGE)

        # check that the link works
        tc.follow('login')

        # auto login for user no 3 
        tc.go('/test/login/3/')
        tc.code(200)

        # user 3 is Fabio
        tc.find('Fabio')
        
        # there is a logout link
        tc.find('logout')

        # check his user profile
        tc.follow ('Fabio')

__test__ = { "doctest": """
>>> c = Client()
"""}

