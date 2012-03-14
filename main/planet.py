"""
Planet Feed collection and 
"""
import os, sys, datetime, urllib
    
from django.conf import settings
from main.server import models
from main.server.const import *

DEBUG = 1

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))
    
def add_feed():
    pass

def download():
    blogs = models.Blog.objects.all()
    
    for blog in blogs:
        fname  = path(settings.PLANET_DIR, "feed-%s.xml" % blog.id)
        try:
            print '*** opening %s, %s' % (blog.id, blog.url)
            stream = open(fname, 'wt')
            data  = urllib.urlopen(blog.url).read()
            stream.write(data)
            stream.close()
        except Exception, exc:
            print '(!) error %s' % exc

if __name__ == '__main__':
    import doctest, optparse
    
    if DEBUG:
        #sys.argv.extend( "--download".split() )
        pass
    
    parser = optparse.OptionParser()
    parser.add_option("--add", dest="url", help="feed url to add to the database", default=None)
    parser.add_option("--download", dest="download", help="downloads feeds", action="store_true", default=False)
    parser.add_option("--update", dest="update", help="adds new posts from downloads", type=int, default=0)
   
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not (opts.url or opts.download or opts.update):
        parser.print_help()
        sys.exit()
        
    if opts.download:
        download()