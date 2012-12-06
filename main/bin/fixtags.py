"""
Fixes tags on the live server
"""
import os, sys, datetime, urllib, glob, re
from itertools import *

from django.conf import settings
from main.server import models, html
from main.server.const import *

def fix_lastuser():
    print "*** fixing last user"
    posts = models.Post.objects.filter(type__in=POST_TOPLEVEL, answer_count__gt=0).order_by("-creation_date")[:limit]
    for post in posts:
        print post.id
        answer = models.Post.objects.filter(parent=post, type=POST_ANSWER).order_by("-lastedit_date")[0]
        post.lastedit_user = answer.author
        post.save()

def lowercase_tags():
    print '*** lowercasing tags'
    posts = models.Post.objects.filter(type__in=POST_TOPLEVEL).order_by("lastedit_date")[:limit]
    for post in posts:
        curr = post.tag_val.lower()
        curr = curr.replace(",", " ")
        if curr != post.tag_val:
            print "saving %s" % post.id
            post.tag_val = curr
            post.set_tags()


def delete_rare():
    tags = models.Tag.objects.filter(count__lt=3)

    for tag in tags:

        posts = tag.post_set.all()

        for post in posts:
            print tag.name, tag.count, post.id
            post.tag_val = post.tag_val.replace(tag.name, "")
            post.set_tags()

        tag.delete()


def run(limit=None):
    fix_lastuser()
    lowercase_tags()
    delete_rare()

if __name__ == '__main__':
    limit = None
    run(limit=limit)