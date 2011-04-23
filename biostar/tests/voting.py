"""
Main test script is executed when running::

    biostar.sh test"

"""
import sys, os, django

from django.test import TestCase
from django.utils import unittest
from django.test.client import Client
from django.conf import settings

from biostar.server import models

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

fixture = path(os.path.dirname(__file__), '..', '..', 'home', 'import', 'test-fixture.json' )

# Post and user ids to use for tests
POST_ID = 1
USER_ID = 3


class VotingTest(TestCase):
    "Tests basic voting and moderation views"

    fixtures = [ fixture ]


    def setUp(self):
        " Login as an admin user"
        self.client = Client()
        r = self.client.post('/admin/password/override/', {'uid':USER_ID, 'password':settings.SECRET_KEY, 'set_type':models.USER_ADMIN})

    def test_vote_up(self):
        "Test addition and removal of an up vote"
        old_score = models.Post.objects.get(pk=POST_ID).score
        
        r = self.client.post('/vote/', {'post':POST_ID, 'type':models.VOTE_UP})
        self.assertContains(r, 'Upvote added')
        score = models.Post.objects.get(pk=POST_ID).score
        self.assertEqual(score, old_score+1)
        
        r = self.client.post('/vote/', {'post':POST_ID, 'type':models.VOTE_UP})
        self.assertContains(r, 'Upvote removed')
        score = models.Post.objects.get(pk=POST_ID).score
        self.assertEqual(score, old_score)
        
    def test_vote_exclusion(self):
        "Tests that up and down votes are mutually exclusive"
        
        old_score = models.Post.objects.get(pk=POST_ID).score
        
        r = self.client.post('/vote/', {'post':POST_ID, 'type':models.VOTE_UP})
        self.assertContains(r, 'Upvote added')
        score = models.Post.objects.get(pk=POST_ID).score
        self.assertEqual(score, old_score+1)
        
        r = self.client.post('/vote/', {'post':POST_ID, 'type':models.VOTE_DOWN})
        self.assertContains(r, 'Downvote added')
        score = models.Post.objects.get(pk=POST_ID).score
        self.assertEqual(score, old_score-1) # Should undo the up and add a down vote
        
    def test_mod_close(self):
        "Tests closing and reopening a post"
        
        
        r = self.client.post('/moderate/', {'post':POST_ID, 'action':'close'})
        self.assertContains(r, 'close performed')
        self.assertTrue(models.Post.objects.get(pk=POST_ID).closed)
        
        r = self.client.post('/moderate/', {'post':POST_ID, 'action':'reopen'})
        self.assertContains(r, 'reopen performed')
        self.assertFalse(models.Post.objects.get(pk=POST_ID).closed)
        
    def test_mod_permissions(self):
        "Tests moderation requiring permissions"
        
        
        # Log in as a non-admin
        r = self.client.post('/admin/password/override/', {'uid':3, 'password':settings.SECRET_KEY, 'set_type':models.USER_NORMAL})
        
        r = self.client.post('/moderate/', {'post':POST_ID, 'action':'close'})
        self.assertContains(r, 'do not have permission')
        self.assertFalse(models.Post.objects.get(pk=POST_ID).closed)
        
    

def suite():
    s = unittest.TestLoader().loadTestsFromTestCase(VotingTest)
    return s
