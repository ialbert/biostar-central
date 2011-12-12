"""
Model testing
"""
import sys, os, django, random, hashlib

from django.utils import unittest
from main.server.const import *
from main.server import models
from main.server.models import User, Post, Tag
from django.test.client import Client

POST_SHOW_URL   = "/post/show/%s/"
NEW_ANSWER_URL  = "/new/answer/%s/"
NEW_COMMENT_URL = "/new/comment/%s/"

def make_uuid():
    "Returns a unique id"
    x = random.getrandbits(256)
    u = hashlib.md5(str(x)).hexdigest()
    return u
   
def create_new_question(c, tag_val='abc efg'):
   
   title   = "Title %s" % make_uuid()
   content = 'Question %s' % make_uuid()
   r = c.post("/new/question/", {'title':title , 'content':content , 'tag_val':tag_val, 'type':POST_QUESTION})
   assert r.status_code == 302
   
   # now select this very question
   q1 = models.Post.all_posts.filter(title=title)[0]
   q2 = models.Post.all_posts.filter(type=POST_QUESTION).order_by('-magic')[0]
   assert q1 == q2
   
   return q1

def note_count(user, unread=False):
   return models.Note.objects.filter(target=user, unread=unread).count()
   
def create_new_answer(c, q):
   
   post_show_url  = POST_SHOW_URL % q.id
   new_answer_url = NEW_ANSWER_URL % q.id
   
   r = c.get(post_show_url)
   assert r.status_code == 200
   
   content = 'Content %s' % make_uuid()
   r = c.post(new_answer_url, {'content':content} )
   assert r.status_code == 302
   
   a1 = models.Post.all_posts.filter(content=content)[0]
   ax = list(q.children.order_by('creation_date'))
   
   # this is the last answer
   assert a1 == ax[-1]
   
   return a1

def create_new_comment(c, obj):
   
   post_show_url  = POST_SHOW_URL % obj.id
   new_comment_url = NEW_COMMENT_URL % obj.id
   
   r = c.get(post_show_url)
   assert r.status_code == 200
   
   content = 'Content %s' % make_uuid()
   r = c.post(new_comment_url, {'content':content} )
   assert r.status_code == 302
   

   
class BiostarUser(unittest.TestCase):
   
   def setUp(self):
      
      # get rid of all notifications
      models.Note.objects.all().delete()
      
      user1, flag1 = User.objects.get_or_create(username='john')
      user1.set_password('test')
      user1.save()
      
      self.c1 = Client()
      self.c1.login(username='john', password='test')
      
      user2, flag1 = User.objects.get_or_create(username='jane')
      user2.set_password('test')
      user2.save()
      self.c2 = Client()
      self.c2.login(username='jane', password='test')
      
   def test_post_answer(self):
      true, eq = self.assertTrue, self.assertEqual
      q  = create_new_question(self.c1)
      
      n1 = note_count(q.author, unread=True)
      a1 = create_new_answer(self.c2, q=q)      
      n2 = note_count(q.author, unread=True)
      eq(n1+1, n2) # the question author gets a new note
      
      n1 = note_count(a1.author, unread=False)
      a2 = create_new_answer(self.c2, q=q)
      n2 = note_count(a1.author, unread=False)
      eq(n1+1, n2)  # the answer author gets a note
      
      a3 = create_new_answer(self.c2, q=q)
      
   def test_post_comment(self):
      true, eq = self.assertTrue, self.assertEqual
      q = create_new_question(self.c1)
      c = create_new_comment(self.c2, obj=q)
      
class BiostarSite(unittest.TestCase):
   
   def test_search(self):
         "Test searching"
         true, eq = self.assertTrue, self.assertEqual

         c = Client()
         r = c.post("/", { 'q':'motif' } )
         eq(r.status_code, 200)
         true( 'motif' in r.content )
   
   def test_answers(self):
         "Test searching"
         true, eq = self.assertTrue, self.assertEqual

         c = Client()
         r = c.post("/", { 'q':'motif' } )
         eq(r.status_code, 200)
         true( 'motif' in r.content )
         
   def test_anav(self):
        "Navigation by unregistered user"
        true, eq = self.assertTrue, self.assertEqual

        c = Client()
        
        locs = "/ /about/ /user/list/ /tag/list/ /badge/list/ /question/unanswered/ ".split()
        for loc in locs:
            r = c.get(loc)
            eq(r.status_code, 200)
        
        r = c.get("/new/question/")
        eq(r.status_code, 302)

   def test_user(self):
        "Visiting a user page by a logged in user"
        true, eq = self.assertTrue, self.assertEqual

        
        user1, flag1 = User.objects.get_or_create(username='john')
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
        user1, flag1 = User.objects.get_or_create(first_name='Jane', last_name='Doe', username='jane', email='jane')
        user1.set_password('test')
        user1.save()
  
        r = c.get("/user/list/")
        eq(r.status_code, 200)
        true('Doe' in r.content)

        c.login(username='jane', password='test')
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

        r = c.get("/")
        eq(r.status_code, 200)
        true(title in r.content)
      

def suite():
    
    s1 = unittest.TestLoader().loadTestsFromTestCase(BiostarUser)
    s2 = unittest.TestLoader().loadTestsFromTestCase(BiostarSite)
    
    suite = unittest.TestSuite( [s1, s2])
    
    return suite
