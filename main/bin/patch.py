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
from datetime import datetime, timedelta

def apply_positive():
    "Removes negative votes from the system"
    posts = models.Post.objects.filter(score__lt=0)
    for post in posts:
        votes = models.Vote.objects.filter(post=post, type=VOTE_DOWN)

        print "Removing %s votes for %s" % ( len(votes), post.id)
        for vote in votes:
            vote.delete()

        post.score = 0
        post.save()


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
        
def blog_cleanup():
    "Updates the bookmark counters. Used after migrating to version 1.2.1"

    # remove banned users
    profs = models.UserProfile.objects.filter(status=USER_BANNED).exclude(about_me="banned")
    if profs:
        profs.update(about_me="banned", website="")

    # move posts for banned users into the blog section (we need a trash section)
    posts = models.Post.objects.filter(author__profile__status=USER_BANNED)
    print "setting up %s posts for deletion" % len(posts)
    if posts:
        posts.update(type=POST_BLOG)

    blogs = models.Post.objects.filter(type=POST_BLOG, status=POST_DELETED)
    
    for blog in blogs:
        print blog.title.encode("ascii", errors="replace")
        try:
            blog.delete()
        except Exception, exc:
            print exc



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
    
def reduce_notes(weeks=30, limit=500):
 
    since = datetime.now() - timedelta(weeks=weeks)
    query = models.Note.objects.filter(date__lt=since)
    dsize = query.count()    
    
    print "*** deleting %s entries" % dsize
    query.delete()
    
    
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
    
def traffic_cleanup(days=1):
    past = datetime.now() - timedelta(days=days)
    query = models.PostView.objects.filter(date__lt=past)

    print "*** deleting %s views" % query.count()

    query.delete()

if __name__ == '__main__':
    import doctest, optparse
  
    #sys.argv.append( '--bookmarks' )
    
    # options for the program
    parser = optparse.OptionParser()
    parser.add_option("-n", dest="n", help="database rows to fetch default=%default", type=int, default=1000)
    
    parser.add_option("-s", dest="s", help="skip this many rows default=%default", type=int, default=0)
    
    parser.add_option("--resave_posts", dest="patt", help="resave posts that match pattern", type=str, default="")
    parser.add_option("--reduce_notes", dest="reduce_notes", help="reduce the number of notification to N", action="store_true", default=False)
    parser.add_option("--reapply_ranks", dest="reapply_ranks", help="reapplies ranks to all posts", action="store_true", default=False)
    parser.add_option("--update_domain", dest="update_domain", help="updates the site domain to match the settings", action="store_true", default=False)
    parser.add_option("--bookmarks", dest="bookmarks", help="updates bookmark counts", action="store_true", default=False)
    parser.add_option("--blog_cleanup", dest="blog_cleanup", help="cleans up deleted blogs", action="store_true", default=False)
    parser.add_option("--positive", dest="positive", help="removes negative ratings", action="store_true", default=False)
    parser.add_option("--traffic_cleanup", dest="traffic", help="removes post view entries older than <days>", type=int, default=0)


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
        reduce_notes(weeks=opts.n)
    
    if opts.patt:
        resave_posts(opts.patt, skip=opts.s, limit=opts.n)
    
    if opts.bookmarks:
        update_bookmark_counts()
        
    if opts.blog_cleanup:
        blog_cleanup()

    if opts.positive:
        apply_positive()

    if opts.traffic:
        traffic_cleanup(opts.traffic)

