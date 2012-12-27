"""
Model testing
"""
import sys, os, django

from django.utils import unittest
from main.server.const import *
from main.server import models
from main.server.models import User, Post, Tag

class BiostarModelTest(unittest.TestCase):
   
   def test_post(self):
        true, eq = self.assertTrue, self.assertEqual
        
        # two users
        user1, flag1 = User.objects.get_or_create( username='john')
        user2, flag2 = User.objects.get_or_create( username='jane')
        
        # save a question to the database
        title, content, tag_val = 'My title', 'My content', 'One tWo threE'
        post1 = Post.objects.create(author=user1, type=POST_QUESTION, title=title, content=content, tag_val=tag_val)
        true( (post1.title, post1.content, post1.tag_val) == (title, content, 'one two three') )
    
        # get rid of all tags
        Tag.objects.all().delete()
        
        # no tags yet in the database    
        tags = [ tag.name for tag in Tag.objects.all() ]        
        eq(0, len(tags))
        
        # now adding some tags
        post1.set_tags()
        names = post1.get_tag_names()
        true( names == [ tag.name for tag in Tag.objects.all() ] )
    
        # one post with these tags
        eq(1, Post.objects.filter(tag_set__name=names[0]).count() )
    
        # create a revision
        models.create_revision(post1, author=user1)
        
        post1.title = post1.title + '?'
        post1.lastedit_user = user2
        post1.save()
        
        # two edits so far
        models.create_revision(post1, user1)
        
        post2 = Post.objects.create(author=user2, type=POST_QUESTION, title=title, content=content, tag_val=tag_val)
        post2.set_tags()
        
        # two posts with this tag
        eq(2, Post.objects.filter(tag_set__name=names[0]).count() )
        



def suite():
    
    suite = unittest.TestLoader().loadTestsFromTestCase(BiostarModelTest)
    
    return suite
