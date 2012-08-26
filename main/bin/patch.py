"""
Applies various corrections and patches to the server
"""
import sys

from django.conf import settings
from main.server import models, html
from main.server.const import *
from django.contrib.sites.models import Site
from django.db.models import Avg, Max, Min, Count
from django.db import transaction
from django.db.models import signals
from itertools import *
from collections import defaultdict

def update_bookmark_counts():
    "Updates the bookmark counters. Used after migrating to version 1.2.1"
    votes = models.Vote.objects.filter(type=VOTE_BOOKMARK).select_related('post')
    
    store = defaultdict(int)
    for vote in votes:
        store[vote.post] += 1
    
    for post, value in store.items():
        print post.title, value
        post.book_count = value
        post.save()
        
def update_domain():
    "This is really only needs to be done once per installation"
    site = Site.objects.get(id=settings.SITE_ID)
    print "*** current site domain %s" % site.domain
    if site.domain != settings.SITE_DOMAIN:
        print '--- updating site domain to %s' % settings.SITE_DOMAIN
        site.domain = settings.SITE_DOMAIN
        site.save()
       
def resave_posts(patt, skip=0, limit=1000):
    "This is really only needs to be done once per installation"

    # this is only needed so that we can update in chunks
    # on a low memory server
    lo = skip * limit
    hi = lo + limit
    
    print "*** resaving posts matching pattern '%s' (%s:%s)" % (patt, lo, hi)

    posts = models.Post.objects.filter(html__contains=patt).order_by('-id')[lo:hi]
    
    for post in posts:
        print "resaving %s, %s" % (post.id, post.title)
        post.save()
    
def remove_notes(target, maxcount=1000):
    """Clears the notes  for each user"""
    
    last_valid = models.Note.objects.filter(target=target).order_by('-id').exclude(sender=target)[maxcount]
    clear_rows = models.Note.objects.filter(target=target, id__lt=last_valid.id).exclude(sender=target)
    clear_rows.delete()
    
    new_count = models.Note.objects.filter(target=target).count()
    print '*** cleared notes for user %s to %s' % (target.id, new_count)
    
    
def reduce_notes(maxcount=1000):
    for user in models.User.objects.all():
        note_count = models.Note.objects.filter(target=user).exclude(sender=user).count()
        if note_count > 2 * maxcount:
            remove_notes(user, maxcount=maxcount)

@transaction.commit_manually()    
def reapply_ranks():
    "Applies the new ranking system on all posts"
    
    
    total = models.Post.objects.filter(type__in=POST_TOPLEVEL).count()
    
    print "*** applying ranks for %s posts " % total
        
    posts = models.Post.objects.filter(type__in=POST_TOPLEVEL).order_by('id')
    
    # disconnect post related signals to speed up update
    #signals.pre_save.disconnect( models.verify_post, sender=models.Post )
    #signals.post_save.disconnect( models.finalize_post, sender=models.Post)
    
    counter = count(1)
    for index, post in izip(counter, posts):
        if (index % 250) == 0:
            perc = 100.0 * index/total 
            print "*** commiting %s, %4.1f%%" % (index, perc) 
            transaction.commit()
        post.rank = html.rank(post)
        post.save()
    transaction.commit()
    

if __name__ == '__main__':
    import doctest, optparse
  
    sys.argv.append( '--update_bookmark_count' )
    
    # options for the program
    parser = optparse.OptionParser()
    parser.add_option("-n", dest="n", help="database rows to fetch default=%default", type=int, default=1000)
    
    parser.add_option("-s", dest="s", help="skip this many rows default=%default", type=int, default=0)
    
    parser.add_option("--resave_posts", dest="patt", help="resave posts that match pattern", type=str, default="")
    parser.add_option("--reduce_notes", dest="reduce_notes", help="reduce the number of notification to N", action="store_true", default=False)
    parser.add_option("--reapply_ranks", dest="reapply_ranks", help="reapplies ranks to all posts", action="store_true", default=False)
    parser.add_option("--update_domain", dest="update_domain", help="updates the site domain to match the settings", action="store_true", default=False)
    parser.add_option("--update_bookmark_counts", dest="update_bookmark_counts", help="updates bookmark counts", action="store_true", default=False)
   
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    if opts.update_domain:
        update_domain()
        
    if opts.reapply_ranks:
        reapply_ranks()
        
    if opts.reduce_notes:
        reduce_notes(maxcount=opts.n)
    
    if opts.patt:
        resave_posts(opts.patt, skip=opts.s, limit=opts.n)
    
    if opts.update_bookmark_counts:
        update_bookmark_counts()