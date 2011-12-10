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
   
   def test_site(self):
        true, eq = self.assertTrue, self.assertEqual
        c = Client()
        
        locs = "/ /about/ /user/list/ /tag/list/ /badge/list/ /question/unanswered/".split()
        for loc in locs:
            r = c.get(loc)
            eq(r.status_code, 200)
        
        r = c.get("/new/question/")
        eq(r.status_code, 302)

        # make a new user
        user1, flag1 = User.objects.get_or_create(first_name='John', last_name='Doe', username='john', email='john')
        user1.set_password('test')
        user1.save()
        
        c.login(username='john', password='test')
        r = c.get("/new/question/")
        eq(r.status_code, 200)
        
        # post a question
        r = c.post("/new/question/", {'title': 'My new question', 'content': 'this is the content of my question', 'tag_val':'aaa bbb ccc'})
        eq(r.status_code, 200)

def suite():
    
    suite = unittest.TestLoader().loadTestsFromTestCase(BiostarSiteTest)
    
    return suite
