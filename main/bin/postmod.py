"""
Modifies an exsisting post in the database
"""
import os, sys, datetime, urllib, glob, string, pprint
    
from django.conf import settings
from main.server import models, html
from main.server.const import *

if __name__ == '__main__':
       
    import optparse
    usage = "usage: %prog [options] file1 file2 ..."
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-p", "--pid", dest="pid", type=int, help="The ID of the post that needs to be changed", default=0)
    parser.add_option("--sticky", dest="sticky", help="set the sticky level 0, 1, 2 etc", type=str, default='')
    parser.add_option("--url", dest="url", help="adss a url to the post", type=str, default='')
   
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not opts.pid:
        parser.print_help()
        sys.exit()
   
    post = models.Post.objects.get(pk=opts.pid) 
    print '*** modifying %s (%s): %s' % (post.get_type_display(), post.id, post.title)
    if opts.sticky:
        post.sticky = int(opts.sticky)
        print '*** setting stickiness to %s' % post.sticky
        post.save()
    
    if opts.url:
        post.url = opts.url
        print '*** setting url to %s' % post.url
        post.save()
