"""
Updating full scores 
"""
import os, sys, datetime, urllib, glob, string, pprint
    
from django.conf import settings
from main.server import models, html
from main.server.const import *

def user_add(username, email, name):
    user, flag = models.User.objects.get_or_create(username=username, email=email)
    if flag:
        print '*** creating user %s' % user
        user.profile.type = USER_MEMBER
        user.profile.display_name = name
        user.profile.save()
    else:
        print '*** found user %s' % user

if __name__ == '__main__':    
    import optparse
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-u", "--username", dest="username", help="The internal username")
    parser.add_option("-e", "--email", dest="email", help="The email for the user")
    parser.add_option("-n", "--name", dest="name",  help="The display name for the user")


    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not (opts.username and opts.email and opts.name):
        parser.print_help()
        sys.exit()
        
    user_add(username=opts.username, email=opts.email, name=opts.name)