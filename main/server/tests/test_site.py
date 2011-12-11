"""
Model testing
"""
import sys, os, django

from django.utils import unittest
from main.server.const import *
from main.server import models
from main.server.models import User, Post, Tag
from django.test.client import Client

class BiostarSiteTest(unittest.TestCase):
   
   def test_anav(self):
        "Navigation by unregistered user"
        true, eq = self.assertTrue, self.assertEqual

        c = Client()
        
        

        locs = "/about/ /user/list/ /tag/list/ /badge/list/ /question/unanswered/".split()
        for loc in locs:
            r = c.get(loc)
            eq(r.status_code, 200)
        
        r = c.get("/new/question/")
        eq(r.status_code, 302)

   def test_user(self):
        "Visiting a user page by a logged in user"
        true, eq = self.assertTrue, self.assertEqual

        
        user1, flag1 = User.objects.get_or_create(first_name='John', last_name='Doe', username='john', email='john')
        user1.set_password('test')
        user1.save()
        
        c = Client()
        c.login(username='john', password='test')
        r = c.get("/user/show/%s/" % user1.id)
        eq(r.status_code, 200)
        r = c.get("/user/edit/%s/" % user1.id)
        eq(r.status_code, 200)

        r = c.post("/user/edit/%s/" % user1.id, { 'display_name':'John', 'email':'abc', 'about_me':'ok', 'location':'xyz', 'website':'xyz'} )
        eq(r.status_code, 302)
        user = User.objects.get(username='john')
        eq(user.email, 'abc')
        eq(user.profile.location, 'xyz')
        eq(user.profile.about_me, 'ok')
      
        r = c.get("/post/list/%s/" % user1.id)
        eq(r.status_code, 200)

        r = c.get("/post/list/%s/questions/" % user1.id)
        eq(r.status_code, 200)

        r = c.get("/post/list/%s/answers/" % user1.id)
        eq(r.status_code, 200)

        r = c.get("/post/list/%s/comments/" % user1.id)
        eq(r.status_code, 200)

   def test_site(self):
        true, eq = self.assertTrue, self.assertEqual
        c = Client()
        
        # make a new user
        user1, flag1 = User.objects.get_or_create(first_name='John', last_name='Doe', username='john', email='john')
        user1.set_password('test')
        user1.save()
  
        r = c.get("/user/list/")
        eq(r.status_code, 200)
        true('Doe' in r.content)

        c.login(username='john', password='test')
        r = c.get("/new/question/")
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
        r = c.post("/new/question/", {'title':title , 'content':content , 'tag_val':'', 'type':POST_QUESTION})
        eq(r.status_code, 200)


        r = c.post("/new/question/", {'title':title , 'content':content , 'tag_val':'aaa bbb ccc', 'type':POST_QUESTION})
        eq(r.status_code, 302)

        r = c.post("/new/question/", {'title':title , 'content':content , 'tag_val':'aaa bbb ccc', 'type':POST_QUESTION}, follow=True)
        eq(r.status_code, 200)

        r = c.post("/preview/", {'content':content})
        eq(r.status_code, 200)

        #r = c.get("/")
        #eq(r.status_code, 200)
        #true(title in r.content)
      

def suite():
    
    suite = unittest.TestLoader().loadTestsFromTestCase(BiostarSiteTest)
    
    return suite
