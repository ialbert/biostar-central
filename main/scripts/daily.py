import sys

from django.conf import settings
from main.server import models, html
from main.server.const import *
from django.contrib.sites.models import Site

from django.db.models import Avg, Max, Min, Count

def update_domain():
    "This is really only needs to be done once per installation"
    site = Site.objects.get(id=settings.SITE_ID)
    print "*** current site domain %s" % site.domain
    if site.domain != settings.SITE_DOMAIN:
        print '--- updating site domain to %s' % settings.SITE_DOMAIN
        site.domain = settings.SITE_DOMAIN
        site.save()
        
def remove_notes(target, maxcount=1000):
    """Clears the notes  for each user"""
    
    last_valid = models.Note.objects.filter(target=target).order_by('-id').exclude(sender=target)[maxcount]
    clear_rows = models.Note.objects.filter(target=target, id__lt=last_valid.id).exclude(sender=target)
    clear_rows.delete()
    
    new_count = models.Note.objects.filter(target=target).count()
    print '*** cleared notes for user %s to %s' % (target.id, new_count)
    
    
def reduce_notelist(maxcount=1000):
    for user in models.User.objects.all():
        note_count = models.Note.objects.filter(target=user).exclude(sender=user).count()
        if note_count > 2 * maxcount:
            remove_notes(user, maxcount=maxcount)
    
def reapply_rank():
    "Applies the new ranking system on all posts"
    posts = models.Post.objects.all()
    
    for post in posts:
        post.rank = html.rank(post)
        post.save()
        
    print "*** reapplied ranks for %s posts " % len(posts)
    
def run():
    #trim_notelist(maxcount=1000)
    reapply_rank()
    pass
    
if __name__ == '__main__':
    import doctest, optparse
   
    # for debugging
    #sys.argv.extend( ["-p", "se0"] )
    
    # options for the program
    parser = optparse.OptionParser()
    parser.add_option("-n", dest="n", help="limit value default=%default", type=int, default=1000)
    parser.add_option("--reduce", dest="reduce", help="reduce the number of notification to N", action="store_true", default=False)
    parser.add_option("--rank", dest="rank", help="reapplies ranks to all posts", action="store_true", default=False)
    parser.add_option("--update_domain", dest="update_domain", help="updates the site domain", action="store_true", default=False)
   
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    if opts.update_domain:
        update_domain()
        
    if opts.rank:
        reapply_rank()
        
    if opts.reduce:
        reduce_notelist(maxcount=opts.n)
    
