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

def merge_tags(fname):
    """Merges tags based on pairs in a file"""
    for line in file(fname):

        val1, val2 = line.strip().split()

        tag1 = models.Tag.objects.filter(name=val1)
        tag2, flag = models.Tag.objects.get_or_create(name=val2)

        if tag1:
            tag1 = tag1[0]

            print "*** merging %s (%s) with %s (%s)" % (tag1.name, tag1.count, tag2.name, tag2.count)
            posts = tag1.post_set.all()
            for post in posts:
                post.tag_val = post.tag_val.replace(tag1.name, tag2.name)
                post.set_tags()
            tag2 = models.Tag.objects.get(name=val2)
            print '*** after merge %s (%s)' % (tag2.name, tag2.count)
            tag1.delete()
def drop_tags():
    "Removes all tags with no open posts"
    for tag in models.Tag.objects.all():
        posts = tag.post_set.all().count()
        if not posts:
            print "*** removing tag %s " % tag.name
            tag.delete()

def remove_tags(fname):

    for line in file(fname):
        name = line.strip()
        tags = models.Tag.objects.filter(name=name)
        if tags:
            tag = tags[0]
            print "*** deleting %s (%s)" % (tag.name, tag.count)
            posts = tag.post_set.all()
            for post in posts:
                post.tag_val = post.tag_val.replace(tag.name, "")
                post.set_tags()
            tag.delete()

def run(limit=None):
    #fix_lastuser()
    #lowercase_tags()
    #delete_rare()
    #old_value, new_value = sys.argv[1:3]
    #merge_tags(old_value, new_value)
    pass

if __name__ == '__main__':
# options for the program
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--merge", dest="merge", help="filename, reads tag-pairs from a file, replaces first tag with the second", default=None)
    parser.add_option("--delete", dest="delete", help="filename, reads tag names from a file, removes them from the system", default=None)
    parser.add_option("--drop", dest="drop", help="drops tags with no posts", action="store_true", default=False)


    (opts, args) = parser.parse_args()

    if not(opts.delete or opts.merge or opts.drop):
        parser.print_help()
        sys.exit()

    if opts.drop:
        drop_tags()

    # stop execution if no parameters were specified
    if opts.merge:
        merge_tags(opts.merge)

    if opts.delete:
        remove_tags(opts.delete)
