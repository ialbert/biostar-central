"""
Main test script is executed when running::

    biostar.sh test"

"""

# disable some spurious warnings
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import md5, sha

import sys, os, django
from django.test import TestCase
from django.utils import unittest
from django.conf import settings

import twill
from twill import commands as tc
from django.core.servers.basehttp import AdminMediaHandler
from django.core.handlers.wsgi import WSGIHandler

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

fixture = path(os.path.dirname(__file__), '..', '..', 'home', 'import', 'test-fixture.json' )

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
    
    def test_homepage(self):
        "Tests the homepage access"
        # go to the home page
        tc.go(HOME_PAGE)
        tc.code(200)

        # check that the link works
        tc.follow('Questions')
        tc.follow('Unanswered')
        tc.follow('Tags')
        tc.follow('Members')
        tc.follow('Badges')

        tc.follow('login')

        tc.go(HOME_PAGE)
        # follow the tag guidelines
        tc.go('/question/tagged/guidelines/')

       


    def test_user(self):
        "Tests a user interaction with a site"

        # this page is active during testing can be turned
        # off for deplyed sites
        tc.go('/admin/password/override/')
        tc.code(200)
        
        # this is how to show what forms are on a pge
        #tc.showforms()

        # submit the override login
        tc.fv(2, 'uid', '3')
        tc.fv(2, 'password', settings.SECRET_KEY)
        tc.submit('submit')
        tc.code(200)

        # user 3 is Fabio
        tc.find('Fabio')
        
        # there is a logout link
        tc.find('logout')

        # check his user profile
        tc.follow ('Fabio')

        # this is how you can save a page for inspection
        #tc.save_html('fabio.html')
        
       
        # in the test dataset he has a reputation of 41
        tc.find('41')

        # on his profile find the title of his question
        tc.find('Finding common motifs')

        # let's ask a question
        tc.go('/question/new/')

        title, content, tags = "Fabios World!", "What is it?", "fabio word"
        tc.fv(2, 'title', title)
        tc.fv(2, 'content', content)
        tc.fv(2, 'tags', tags)
        tc.submit('submit')
        tc.code(200)

        #tc.save_html('test.html')

        # the resulting page will have the question
        tc.find(title)
        tc.find(content)

        # check that it appeared in the questions list
        tc.follow('Questions')
        tc.find(title)
        tc.follow(title)

        # user edits their question
        tc.follow('edit')
        # edit the question
        title, content, tags = "Fabios WorldXYZ!", "What is itXYZ?", "fabioX wordY"
        tc.fv(2, 'title', title)
        tc.fv(2, 'content', content)
        tc.fv(2, 'tags', tags)
        tc.submit('submit')
        tc.find(title)
        tc.find(content)

        #tc.save_html('test.html')
        #tc.showforms()
        # let's add an answer
       
        # now a second user logs in
        tc.go('/admin/password/override/')
        tc.code(200)
        tc.fv(2, 'uid', '2')
        tc.fv(2, 'password', settings.SECRET_KEY)
        tc.submit('submit')
        tc.code(200)

        # finds Faboios question
        tc.follow('Questions')
        tc.find(title)
        tc.follow(title)

        tc.showforms()

        # answers Fabios quesion
        content = 'Read the manual FABIO!!!!'
        # why is this in twill form no 3
        tc.fv(3, 'content', content)
        tc.submit('submit')
        tc.code(200)
        tc.find(content)

        # also upvotes the question since he answered it ;-)
        tc.go("/vote/%s/preview/" % HOME_PAGE)

        # check that markdown preview page works
        tc.go("%s/preview/" % HOME_PAGE)
        tc.code(200)
        tc.find("no input")

 
        


def suite():
    s = unittest.TestLoader().loadTestsFromTestCase(FunctionalTest)
    return s
