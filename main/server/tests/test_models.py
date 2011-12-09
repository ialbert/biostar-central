"""
Model testing
"""
import sys, os, django

from django.utils import unittest

"""
  >>> post.set_tags()
    
    >>> user2, flag = User.objects.get_or_create(first_name='Jane', last_name='Doe', username='jane', email='jane')
    >>>
    >>> title, content, tag_val = 'My title', 'My content', 'One tWo threE'
    >>> (post1.title, post1.content, post1.tag_val) == (title, content, tag_val)
    True
    >>> [ tag.name for tag in Tag.objects.all() ]
    []
    >>> post1.set_tags()
    >>> names = post1.get_tag_names()
    >>> names == [ tag.name for tag in Tag.objects.all() ]
    True
    >>> Post.objects.filter(tag_set__name=names[0]).count()
    1
    >>> post1.revisions.count()
    1
    >>> post2 = Post.objects.create(author=user, type=POST_QUESTION, title=title, content=content, tag_val=tag_val)
    >>> post2.set_tags()
    >>> Post.objects.filter(tag_set__name=names[0]).count()
    2
    >>> post2.title = post2.title + '?'
    >>> post2.lastedit_user = user2
    >>> post2.save()
    >>> post2.revisions.count()
    2
    
    
    
   
    True
    
    >>> post1.set_tags()
    >>> names = post1.get_tag_names()
    >>> names == [ tag.name for tag in Tag.objects.all() ]
    True
    
    
    >>> post2 = Post.objects.create(author=user, type=POST_QUESTION, title=title, content=content, tag_val=tag_val)
    >>> post2.set_tags()
    >>> Post.objects.filter(tag_set__name=names[0]).count()
    2
    
    >>> post2.title = post2.title + '?'
    >>> post2.lastedit_user = user2
    >>> post2.save()
    >>> post2.revisions.count()
    2
"""

from main.server.const import *
from main.server.models import User, Post, Tag

class BiostarModelTest(unittest.TestCase):
   
   def test_post(self):
        true, eq = self.assertTrue, self.assertEqual
        
        # two users
        user1, flag1 = User.objects.get_or_create(first_name='John', last_name='Doe', username='john', email='john')
        user2, flag2 = User.objects.get_or_create(first_name='Jane', last_name='Doe', username='jane', email='jane')
        
        # save a question to the database
        title, content, tag_val = 'My title', 'My content', 'One tWo threE'
        post1 = Post.objects.create(author=user1, type=POST_QUESTION, title=title, content=content, tag_val=tag_val)
        true( (post1.title, post1.content, post1.tag_val) == (title, content, tag_val) )
    
        # no tags yet in the database    
        tags = [ tag.name for tag in Tag.objects.all() ]        
        eq(0, len(tags))
        
        # now adding some tags
        post1.set_tags()
        names = post1.get_tag_names()
        true( names == [ tag.name for tag in Tag.objects.all() ] )
    
        # one post with these tags
        eq(1, Post.objects.filter(tag_set__name=names[0]).count() )
    
        # one revision so far
        true(1 == post1.revisions.count())
        post1.title = post1.title + '?'
        post1.lastedit_user = user2
        post1.save()
        
        # two edits so far
        eq(2, post1.revisions.count())
        post2 = Post.objects.create(author=user2, type=POST_QUESTION, title=title, content=content, tag_val=tag_val)
        post2.set_tags()
        
        # two posts with this tag
        eq(2, Post.objects.filter(tag_set__name=names[0]).count() )
        
        



def suite():
    
    suite = unittest.TestLoader().loadTestsFromTestCase(BiostarModelTest)
    
    return suite
