"""
Planet Feed collection and 
"""
import os, sys, datetime, urllib, glob
    
from django.conf import settings
from main.server import models, html
from main.server.const import *
import feedparser

DEBUG = False

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

def get_datapath(blog):
    return path(settings.PLANET_DIR, "feed-%s.xml" % blog.id)

def init(limit):
    "Initialize blogs to a few known blogs"
    
    urls =[
        'http://plindenbaum.blogspot.com/feeds/posts/default',
        'http://nsaunders.wordpress.com/feed/',
        'http://feeds2.feedburner.com/bcbio',
        'http://ivory.idyll.org/blog/tags/bioinformatics?flav=atom',
        'http://feeds.feedburner.com/GenomesUnzipped?format=xml',
        'http://feeds.feedburner.com/Massgenomics?format=xml',
        'http://feeds.feedburner.com/MyWeblogOnBioinformaticsGenomeScienceNextGenerationSequencing',
        'http://feeds.feedburner.com/OmicsOmics',
        'http://feeds.feedburner.com/JermdemoRaisedToTheLaw?format=xml',
        'http://hackmap.blogspot.com/feeds/posts/default',
        'http://feeds.feedburner.com/GettingGeneticsDone?format=xml',
        'http://genomeinformatician.blogspot.com/feeds/posts/default',
    ]
    
    for url in urls[:limit]:
        add(url)
        
def add(url):
    
    try:
        print '*** adding %s' % url
        fname = path(settings.PLANET_DIR, 'add-url.xml')
        text = urllib.urlopen(url).read()
        stream = file(fname, 'wt')
        stream.write(text)
        stream.close()
        doc   = feedparser.parse(fname)
        title = doc.feed.title
        desc  = doc.feed.description
        
        username = title[:30]
        
        users = models.User.objects.filter(username=username)
        
        if users:
            print '(!) blog name name %s exists' % username
            return
        
        print '*** adding user %s' % username
        user = models.User.objects.create_user(username=username, email=username)
        user.profile.display_name = title
        user.profile.type = USER_BLOG
        user.profile.website  = url
        user.profile.about_me = desc
        user.profile.save()
        blog = models.Blog(author=user, url=url)
        blog.save()
        
    except Exception, exc:
        print '(!) error %s' % exc

def title(row):
    title = row.title[:200]
    return title

def update(limit):
    blogs = models.Blog.objects.all()
    for blog in blogs:
        posts = models.Post.objects.filter(author=blog.author)
        seen  = set( p.title for p in posts )
        fname = get_datapath(blog)
        try:
            print '*** reading %s, %s' % (blog.id, blog.url)
            now = datetime.datetime.now()
            models.UserProfile.objects.filter(user=blog.author).update(last_visited=now)
            doc = feedparser.parse(fname)
            
            ent = [ e for e in doc.entries if title(e) not in seen ]
            ent = ent[:limit]
            for r in ent:
                if not r.title:
                    continue;
                date = r.date_parsed
                date = datetime.datetime(date[0], date[1], date[2])
                content = html.strip_tags(r.description)
                post = models.Post(title=title(r), url=r.link, author=blog.author,  type=POST_BLOG, content=content, creation_date=date)
                post.save()
                print '*** added post %s' % post.title.encode("ascii", errors='replace')
        except KeyError, exc:
            print '(!) error %s' % exc
    
def download():
    blogs = models.Blog.objects.all()    
    for blog in blogs:
        fname  = get_datapath(blog)
        try:
            print '*** downloading blog %s at %s' % (blog.id, blog.url)
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
        #sys.argv.extend( "--update 2".split() )
        sys.argv.extend( "--url http://plindenbaum.blogspot.com/feeds/posts/default".split() )
    
    parser = optparse.OptionParser()
    parser.add_option("--url", dest="url", help="feed url to add to the database", default=None)
    parser.add_option("--download", dest="download", help="downloads feeds", action="store_true", default=False)
    parser.add_option("--update", dest="update", help="adds new posts from downloads", type=int, default=0)
    parser.add_option("--init", dest="init", help="initializes feeds", type=int, default=0)
    
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not (opts.url or opts.download or opts.update or opts.init):
        parser.print_help()
        sys.exit()
    
    if opts.init:
        init(opts.init)
        
    if opts.url:
        add(opts.url)
        
    if opts.download:
        download()
        
    if opts.update:
        update(opts.update)