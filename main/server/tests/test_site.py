"""
Model testing
"""
import sys, os, django, random, hashlib
from django.conf.urls import url, patterns, include
from django.core.urlresolvers import reverse
from django.conf import settings

from django.utils import unittest
from main.server.const import *
from main.server import models, search
from main.server.models import User, Post, Tag
from django.test.client import Client
from django.test import TestCase

def user_show(c, uid, post_type=None):
    if post_type:
        loc = reverse("show-user-content", kwargs={'uid':uid, 'post_type':post_type})
    else:
        loc = reverse("show-user", kwargs={'uid':uid})
    r = c.get(loc)
    return r

def user_profile(c, uid, tab=None):
    if tab:
        loc = reverse("user-profile-tab", kwargs={'uid':uid, 'tab':tab})
    else:
        loc = reverse("user-profile", kwargs={'uid':uid})
    r = c.get(loc)
    return r


class UserActions(TestCase):
        
    def test_user_edit(self):
        true, eq = self.assertTrue, self.assertEqual
        
        user1, flag1 = User.objects.get_or_create(username='john')
        user1.set_password('test')
        user1.save()
        
        c = Client()
        c.login(username='john', password='test')
        
        uid1 = user1.id
        
        r = user_show(c, uid=uid1)
        eq(r.status_code, 200)
        
        # edit the user
        url = reverse('user-edit', kwargs={'uid':uid1})

        r = c.get(url)
        eq(r.status_code, 200)
        
        r = c.post(url, { 'display_name':'John', 'email':'abc', 'about_me':'ok', 'location':'xyz', 'website':'xyz'} )
        eq(r.status_code, 302)
        
        user = User.objects.get(username='john')
        eq(user.email, 'abc')
        eq(user.profile.location, 'xyz')
        eq(user.profile.about_me, 'ok')

    def test_user_post(self):
        true, eq = self.assertTrue, self.assertEqual
        
        user1, flag1 = User.objects.get_or_create(username='john')
        user1.set_password('test')
        user1.save()
        
        c = Client()
        
        url = reverse("new-post")
        r = c.get(url)
        eq(r.status_code, 302)
        
        c.login(username='john', password='test')

        url = reverse("new-post")
        r = c.get(url)
        eq(r.status_code, 200)
        
        # post a question
        title   = 'My new question about XYZBC123'
        content = '''
        This is the content of my question!
        
            print 123
            print 123

        some other *markup in it*

            for x in range(10):
                print x
        '''

        # missing tag will make it stay on the same page with an error message
        r = c.post(url, {'title':title , 'content':content , 'tag_val':'', 'type':POST_QUESTION})
        eq(r.status_code, 200)
        true('required' in r.content)
        
        # filling in the tags make it move to redirect
        r = c.post(url, {'title':title , 'content':content , 'tag_val':'aaa bbb ccc', 'type':POST_QUESTION})
        eq(r.status_code, 302)
        
        r = c.get("/")
        eq(r.status_code, 200)
        true(title in r.content)
        
class DataNav(TestCase):
    "Navigates to each static page that contains data"
    fixtures = [ 'simple-test.json' ]
    
    def setUp(self):
        search.VERBOSE = 0
        search.full_index()
    
    def test_voting(self):
        true, eq = self.assertTrue, self.assertEqual
        
        user1 = User.objects.get(id=2)
        user1.set_password('test')
        user1.save()
        
        c = Client()
        url = reverse("vote")
    
        # user not logged in
        r = c.post(url, {'post':2 , 'type':'upvote'})
        true("error" in r.content)
        eq(r.status_code, 200)
        
        # voting on self post
        c.login(username=user1.username, password='test')
        r = c.post(url, {'post':2 , 'type':'upvote'})
        true("error" in r.content)
        eq(r.status_code, 200)

        c.login(username=user1.username, password='test')
        r = c.post(url, {'post':2 , 'type':'bookmark'})
        true("success" in r.content)
        eq(r.status_code, 200)


        # correct voting
        r = c.post(url, {'post':4 , 'type':'upvote'})
        true("success" in r.content)
        eq(r.status_code, 200)
        
    def test_uid(self):
        "Navigation to pages that require UID parameters"
        true, eq = self.assertTrue, self.assertEqual

        c = Client()
        
        # show all activity for a user
        r = user_show(c, uid=2)
        eq(r.status_code, 200)
        
        # test user post content
        r = user_show(c, uid=2, post_type="answer")
        eq(r.status_code, 200)
        true("SHRiMP" in r.content)
        
        # test user profile
        r = user_profile(c, uid=2)
        eq(r.status_code, 200)
        
        # test bookmarks
        r = user_profile(c, uid=2, tab="bookmarks")
        eq(r.status_code, 200)
        true("no bookmarks" in r.content)
        
        # test badges
        r = user_profile(c, uid=2, tab="badges")
        eq(r.status_code, 200)
        true("Supporter" in r.content)
        
        # test moderator tab
        r = user_profile(c, uid=2, tab="moderator")
        eq(r.status_code, 200)
        true("Deleted" in r.content)
        
        # test searching from the main page
        true, eq = self.assertTrue, self.assertEqual
        r = c.get(reverse("search"), { 'q':'motif' } )
        eq(r.status_code, 200)
        true( 'common motifs' in r.content )
        true( 'Comment' in r.content )
        
        r = c.get(reverse("search"), { 'q':'motif', 't':"Question" } )
        eq(r.status_code, 200)
        true( 'common motifs' in r.content )
        true( 'Comment' not in r.content )
            
class SimpleNav(TestCase):
    "Navigates to each static page"
    
    def test_simple(self):
        "Navigation to standalone pages"
        true, eq = self.assertTrue, self.assertEqual
        c = Client()
        locs = "index about rss faq badge-list tag-list user-list".split()
        for loc in locs:
            r = c.get(reverse(loc))
            eq(r.status_code, 200)
        
        # get a redirect when trying to link to a question
        r = c.get(reverse("new-post"))
        eq(r.status_code, 302)
        
    def test_content(self):
        "Navigation to content pages"
        true, eq = self.assertTrue, self.assertEqual

        c = Client()
        args = "show,questions show,forum show,tutorial show,recent".split()
        for pair in args:
            name, arg = pair.split(",")
            loc = reverse(name, kwargs={'tab':arg})
            r = c.get(loc)
            eq(r.status_code, 200)
            
    def test_create(self):
        true, eq = self.assertTrue, self.assertEqual

        c = Client()
        pass
    
def suite():
    
    simple = unittest.TestLoader().loadTestsFromTestCase(SimpleNav)
    data = unittest.TestLoader().loadTestsFromTestCase(DataNav)
    uact = unittest.TestLoader().loadTestsFromTestCase(UserActions)
    
    suite = unittest.TestSuite( [ uact, simple, data ])
    
    return suite
