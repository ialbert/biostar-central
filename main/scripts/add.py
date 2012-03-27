"""
Updating full scores 
"""
import os, sys, datetime, urllib, glob, string
    
from django.conf import settings
from main.server import models, html
from main.server.const import *

DEBUG = 0
    
def parse(fname):
    """
    An importable post is self contained file with title and tags in it:
    
    TITLE: Post Title Goes Here
    TAGS: tag1 tag2
    The rest is the content of the post
    """
    lines = file(fname).readlines()
    title = lines[0].split(':')[-1]
    tags  = lines[1].split(':')[-1]
    body  = ''.join(lines[2:])
    return map(string.strip, (title, tags, body))
    
def add_files(fnames, uid, ptype):
    user = models.User.objects.get(id=uid)
    
    for fname in fnames:
        title, tag_val, body = parse(fname)
        print '*** adding %s' % title
        post = models.Post(title=title, author=user,  type=ptype, tag_val=tag_val, content=body)
        post.save()
        post.set_tags()

def get_user(name, email, website='', about_me=''):
    users = models.User.objects.filter(email=email)
    if users:
        print '*** email %s already exists' % email
        return users[0]
    else:
        print '*** creating user %s:%s' % (name, email)
        username = models.make_uuid()[:10]
        user = models.User.objects.create(username=username, first_name=name, email=email)
        user.profile.display_name = name    
        user.profile.website  = website
        user.profile.about_me = about_me
        user.profile.save()
        return user

if __name__ == '__main__':
    # debug options for the program
    if DEBUG:
        sys.argv.extend( "-u 2 test/post.txt test/post.txt".split() )
        
    import optparse
    usage = "usage: %prog [options] file1 file2 ..."
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-u", "--userid", dest="uid", help="The ID of the user that generated the post", default=None)
    parser.add_option("-t", "--type", dest="ptype", help="The post type", default='Question')
    
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not opts.uid or not args:
        parser.print_help()
        sys.exit()
        
    ptype = opts.ptype.lower()
    ptype = POST_REV_MAP[ptype]
    add_files(fnames=args, uid=opts.uid, ptype=ptype)