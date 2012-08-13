"""
Updating full scores 
"""
import os, sys, datetime, urllib, glob, string, pprint
    
from django.conf import settings
from main.server import models, html
from main.server.const import *

def usermod(uid, role):
    user = models.User.objects.get(pk=uid)
    user.profile.type = role
    user.profile.save()
    
if __name__ == '__main__':    
    import optparse
    usage = "usage: %prog [options] file1 file2 ..."
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-u", "--userid", dest="uid", type=int, help="The ID of the user", default=0)
    parser.add_option("-r", "--role", dest="role", type=int, help="The role for this user", default=0)
    
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not (opts.uid and opts.role):
        parser.print_help()
        print "*** valid roles"
        pprint.pprint(USER_TYPES)
        sys.exit()
        
    usermod(uid=opts.uid, role=opts.role)