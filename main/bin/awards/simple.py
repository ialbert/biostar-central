"""
Simple (single query awards) 
"""
import os, sys, datetime, urllib, glob, string
    
from django.conf import settings
from main.server import models, html, notegen
from main.server.const import *
from django.db import connection, transaction
from django.db.models import Count

DEBUG = 1
  
def instant(request, user):
    "Produces an instant award if applicable"
    awards = set( models.Award.objects.values('id') )
    print awards
    
def notify(user, award):
    text = notegen.awardnote(award.badge)
    note = models.Note.objects.create(sender=user, target=user, content=text)
    
@transaction.commit_on_success
def apply(badge, users):
    for user in users:
        print '*** awarding %s to %s, %s' % (badge.name, user.id, user.email)
        award = models.Award.objects.create(user=user, badge=badge)
        notify(user=user, award=award)

def teacher(badge):    
    # posted an upvoted answer
    users = models.User.objects.filter(post__type=POST_ANSWER, post__score__gt=0).exclude(award__badge=badge).distinct()
    apply(badge, users)

def supporter(badge):    
    # at least one upvote
    users = models.User.objects.filter(vote__type=VOTE_UP ).exclude(award__badge=badge).distinct()
    apply(badge, users)

def civic_duty(badge, val=300):    
    # users that number of votes
    users = models.User.objects.annotate(vcount=Count('vote')).filter(vcount__gt=val).exclude(award__badge=badge).distinct()
    apply(badge, users)

def stellar(badge, val=0):    
    # question bookmarked by 20 users
    votes = models.Vote.objects.all().annotate(pcount=Count('post'))

    
AWARDS = [
    #('Teacher', teacher),
    #('Supporter', supporter),
    #('Civic Duty', civic_duty),
    ('Stellar Question', stellar),

]

def list_badges():
   badges = get_badges()
   for badge in badges:
        print "%s %s: %s" % (badge.id, badge.name, badge.description)

def award():
    
    #models.Award.objects.all().delete()
    
    for name, func in AWARDS:
        badge  = models.Badge.objects.get(name=name)
        print '*** checking badge %s: %s' % (badge.name, badge.description)
        func(badge)
        
    
if __name__ == '__main__':
    # debug options for the program
    if DEBUG:
        sys.argv.extend( "1 2 3".split() )
        
    import optparse
    usage = "usage: %prog [options] id1 id2 id3 ..."
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--list", dest="list", default=False, action='store_true', help="Lists the awards")
    
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    
    if not (args or opts.list):
        parser.print_help()
        sys.exit()
    
    if opts.list:
        list_badges()
    
    if args:
        award()
    #add(fnames=args, uid=opts.uid, ptype=ptype)