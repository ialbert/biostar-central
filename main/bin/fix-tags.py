"""
Fixes tags on the live server
"""
import os, sys, datetime, urllib, glob, re
from itertools import *

from django.conf import settings
from main.server import models, html
from main.server.const import *

def fix_mytags(text):
    out = []
    text = text.replace(" ", "+")
    text = text.replace(",", "+")
    text = text.replace("++", "+")
    text = text.replace("++", "+")
    text = text.lower()
    return text
  
def needs_fixing(post):
    val = post.tag_val
    con = False
    for char in '-,."(':
        con = con or char in val
    return con

def run(limit=None):
    
    print '--- fixing mytags'
    profs = models.UserProfile.objects.all().select_related('user').all()[:limit]
    for p in profs:
        before = p.my_tags
        p.my_tags = fix_mytags(p.my_tags)
        after  = p.my_tags
        if before != after:
            print p.user, p.user.email, p.my_tags
            p.save()
            print '---'
    
    print '--- fixing posts'
    
    posts = models.Post.objects.filter(type__in=POST_TOPLEVEL).all()[:limit]
    posts = ifilter(needs_fixing, posts)
    for post in posts:
        print "Post %s, %s" %  (post.id, post.tag_val)
        post.set_tags()
        print '---'
    
if __name__ == '__main__':
    run()