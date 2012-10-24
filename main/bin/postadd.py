"""
Imports a new post and adds it to the database 
"""
import os, sys, datetime, urllib, glob, string, pprint
    
from django.conf import settings
from main.server import models, html
from main.server.const import *

uni = html.unicode_or_bust

DEBUG = 0

def error():
    print """
The file to be added as a post must have the form:

TITLE: Post title
TAGS: tag1 tag2 tag3
<post text goes here>
"""
    sys.exit('*** stopping')

def parse(fname):
    """
    An importable post is self contained file with title and tags in it:
    
    TITLE: Post Title Goes Here
    TAGS: tag1 tag2
    The rest is the content of the post
    """
    lines = file(fname).readlines()
    title = lines[0]
    if not title.startswith("TITLE:"):
        error()
    title = title.split(':')[-1]

    tags  = lines[1]
    if not tags.startswith("TAGS:"):
        error()
    tags  = tags.split(':')[-1]

    body = ''.join(lines[2:])
    vals = map(string.strip, (title, tags, body))
    vals = map(uni, vals)

    return vals

def add_files(fnames, uid, ptype, delete=False, sticky=False):
    
    user = models.User.objects.get(pk=uid)
    
    for fname in fnames:
        title, tag_val, body = parse(fname)
        post = models.Post(title=title, author=user,  type=ptype, tag_val=tag_val, content=body, sticky=sticky)
        post.save()
        print '*** adding %s, %s, %s, %s' % (post.id, user.id, user.email, title)
        post.set_tags()
        if delete:
            os.delete(fname)

if __name__ == '__main__':
    # debug options for the program
    if DEBUG:
        sys.argv.extend( "-u 2 test/post.txt test/post.txt".split() )
        
    import optparse
    usage = "usage: %prog [options] file1 file2 ..."
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-u", "--userid", dest="uid", type=int, help="The ID of the user that generated the post", default=2)
    parser.add_option("-e", "--email", dest="email", help="The email of the user that generated the post. Overrides the uid parameter")
    parser.add_option("-t", "--type", dest="ptype", type=int, help="The type of the post", default=1)
    parser.add_option("--delete", dest="delete", help="deletes the imported file", action="store_true", default=False)
    parser.add_option("--sticky", dest="sticky", help="set the sticky flag", action="store_true", default=False)
   
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not (opts.uid or opts.email) or not args:
        parser.print_help()
        pprint.pprint(POST_TYPES)
        sys.exit()

    if opts.email:
        uid = models.User.objects.get(email=opts.email).id
    else:
        uid = opts.uid

    #ptype = POST_REV_MAP[ptype]
    add_files(fnames=args, uid=uid, ptype=opts.ptype, delete=opts.delete, sticky=opts.sticky)