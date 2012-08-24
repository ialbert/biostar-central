"""
Modifies an exsisting post in the database
"""
import os, sys, datetime, urllib, glob, string, pprint
    
from django.conf import settings
from main.server import models, html
from main.server.const import *

def post_mod(pid, **kwds):
    
    post = models.Post.objects.get(pk=pid)
    for key, value in kwds.items():
        setattr(post, key, value)
        print post, post.sticky
    
    post.save()

if __name__ == '__main__':
       
    import optparse
    usage = "usage: %prog [options] file1 file2 ..."
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-p", "--pid", dest="pid", type=int, help="The ID of the post that needs to be changed", default=0)
    parser.add_option("--sticky", dest="sticky", help="set the sticky flag", action="store_true", default=False)
   
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not opts.pid:
        parser.print_help()
        sys.exit()
        
    #ptype = POST_REV_MAP[ptype]
    post_mod(pid=opts.pid, sticky=opts.sticky)