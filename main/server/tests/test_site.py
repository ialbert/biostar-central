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
    fixtures = [ 'simple-test.json.gz' ]
    
    def setUp(self):
        search.VERBOSE = 0
        search.full_index()
    
    def test_post_feeds(self):
        "Tests that all the feeds work"
        true, eq = self.assertTrue, self.assertEqual
        
        c = Client()
        url = reverse("latest-feed")
        r = c.get(url)
        eq(r.status_code, 200)
    
        url = reverse("tags-feed", kwargs={'text':'bwa+blast'})
        r = c.get(url)
        eq(r.status_code, 200)
        
        url = reverse("post-feed", kwargs={'text':'2+4'})
        r = c.get(url)
        eq(r.status_code, 200)
        
        url = reverse("user-feed", kwargs={'text':'2+3+4'})
        r = c.get(url)
        eq(r.status_code, 200)
        
        uuid = models.User.objects.get(id=2).profile.uuid
        url = reverse("mytags-feed", kwargs={'uuid':uuid})
        r = c.get(url)
        eq(r.status_code, 200)
        
        url = reverse("notification-feed", kwargs={'uuid':uuid})
        r = c.get(url)
        eq(r.status_code, 200)
        
     
    
    def test_post_content(self):

        true, eq = self.assertTrue, self.assertEqual

        c = Client()
        url = reverse("post-show", kwargs={'pid':13})
        r = c.get(url)
        eq(r.status_code, 200)
    
        true( "CHIP DNA" in r.content)
        true("Chang" in r.content)
        true("Zhang" in r.content)
        
        post = models.Post.objects.get(id=13)
        url = reverse("post-edit", kwargs={'pid':13})
        
        # unauthorized user trying to edit
        r = c.post(url, {'title':post.title , 'content':post.content , 'tag_val':'', 'type':POST_QUESTION}, follow=True)
        eq(r.status_code, 200)
        true('OpenID' in r.content)
        
        joe = User.objects.get(id=5)
        joe.set_password('test')
        joe.save()
        c.login(username=joe.username, password='test')
        r = c.post(url, {'title':post.title , 'content':post.content , 'tag_val':'', 'type':POST_QUESTION}, follow=True)
        eq(r.status_code, 200)
        true("may not edit" in r.content)
        
        mod = User.objects.get(id=2)
        mod.set_password('test')
        mod.save()
        c.login(username=mod.username, password='test')
        
        # missing tag, will not submit the post
        title = "ABCDEFG"
        r = c.post(url, {'title':title, 'content':post.content , 'tag_val':'', 'type':POST_QUESTION})
        eq(r.status_code, 200)
        true(models.Post.objects.get(id=13).title == post.title)
    
        # missing tag, will not submit the post
        r = c.post(url, {'title':title , 'content':post.content , 'tag_val':'ABCD', 'type':POST_QUESTION})
        eq(r.status_code, 302)
        true(models.Post.objects.get(id=13).title == title)
        
        # test the request merge pages
        url = reverse("request-merge")
        r = c.get(url)
        eq(r.status_code, 200)
        
        url = reverse("approve-merge",kwargs={'master_id':2, 'remove_id':3})
        r = c.get(url)
        eq(r.status_code, 302)
        
    def test_moderator(self):
        "Testing moderator actions"
        true, eq = self.assertTrue, self.assertEqual
        
        joe = User.objects.get(id=5)
        joe.set_password('test')
        joe.save()
        
        mod = User.objects.get(id=2)
        mod.set_password('test')
        mod.save()
        
        c = Client()
        
        post_close = reverse("post-moderate", kwargs={'pid':4, 'status':"close"})
        post_open  = reverse("post-moderate", kwargs={'pid':4, 'status':"open"})
        user_suspend = reverse("user-moderate", kwargs={'uid':4, 'status':"suspend"})
        user_reinstate = reverse("user-moderate", kwargs={'uid':4, 'status':"reinstate"})
        
        # trying to moderate as a non moderator
        r = c.get(post_close)
        r = c.get(user_suspend)
        
        true(Post.objects.get(id=4).status == POST_OPEN)
        true(User.objects.get(id=4).profile.status == USER_ACTIVE)
        
        # moderation by a non moderator
        c.login(username=joe.username, password='test')
        r = c.get(post_close)
        r = c.get(user_suspend)
        true(Post.objects.get(id=4).status == POST_OPEN)
        true(User.objects.get(id=4).profile.status == USER_ACTIVE)
        
        # moderation by a moderator
        c.login(username=mod.username, password='test')
        r = c.get(post_close)
        r = c.get(user_suspend)
        true(Post.objects.get(id=4).status == POST_CLOSED)
        true(User.objects.get(id=4).profile.status == USER_SUSPENDED)
        
        # reopen, reinstate users
        r = c.get(post_open)
        r = c.get(user_reinstate)
        true(Post.objects.get(id=4).status == POST_OPEN)
        true(User.objects.get(id=4).profile.status == USER_ACTIVE)
            
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
        
        r = c.get(reverse("search"), { 'q':'motif', 't':"Question" } )
        eq(r.status_code, 200)
        true( 'common motifs' in r.content )
            
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
        tabs = "recent popular mytags questions unanswered planet forum tutorials".split()
        for tab in tabs:
            r = c.get(reverse("show", kwargs={'tab':tab}))
            eq(r.status_code, 200)
            
    def Xtest_post_edit(self):
        true, eq = self.assertTrue, self.assertEqual

        c = Client()
        url = reverse("post-edit", kwargs={'pid':13})
        r = c.get(url)
        eq(r.status_code, 200)

        print r.content
        
        true( "CHIP DNA" in r.content)
        true("Chang" in r.content)
        true("Zhang" in r.content)
        
        
def suite():
    
    simple = unittest.TestLoader().loadTestsFromTestCase(SimpleNav)
    data = unittest.TestLoader().loadTestsFromTestCase(DataNav)
    uact = unittest.TestLoader().loadTestsFromTestCase(UserActions)
    
    suite = unittest.TestSuite( [ uact, simple, data ])
    
    return suite
