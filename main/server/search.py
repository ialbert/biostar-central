"""
Indexes all post content
"""    
    
import shutil, os, gc
from django import forms
from django.conf import settings
from main.server.const import *
from itertools import *
from main.server import html, const, formdef
from main.server.html import get_page
from whoosh import store, fields, index, highlight
from whoosh.qparser import QueryParser,  MultifieldParser, WildcardPlugin
from whoosh.analysis import StemmingAnalyzer, StopFilter, StandardAnalyzer
from django.contrib import messages
from itertools import *

# common words that should be ignored
from main.server.models import Post

stoplist = "read reads lab most post posted posting\
".split()

stem = StemmingAnalyzer(stoplist=stoplist, minsize=3, maxsize=40, cachesize=-1)
stop = StopFilter()
full = stem | stop 

SCHEMA = fields.Schema(
    pid     = fields.NUMERIC(stored=True, unique=True),
    title   = fields.TEXT(analyzer=full, stored=True),
    content = fields.TEXT(analyzer=full, stored=True), 
    type    = fields.TEXT(stored=True),
)

def initialize(sender=None, **kwargs):
    "Initializes the index. Called from a signal"
    if sender.__name__ != 'main.server.models':
        return
    path = settings.WHOOSH_INDEX
    
    if not os.path.exists(path):
        print('*** creating %s' % path)
        os.mkdir(path)
        print ('*** initializing search index %s' % path)
        ix = index.create_in(path, SCHEMA)
    
choices = [
    ("all", "All types"), ("top", "Top level"),  ("Question", "Questions"), ("Tutorial", "Tutorials"), ("Forum", "Forum"),
]

VERBOSE = 1
def info(msg):
    if VERBOSE:
        print "*** %s" % msg
        
class SearchForm(forms.Form):
    "A form representing a new question"
    q = forms.CharField(max_length=200,  initial="", widget=forms.TextInput(attrs={'size':'50', 'class': 'span6', 'placeholder': 'Search Biostar'}))
    t = forms.ChoiceField(choices=choices, required=False)

def safe_int(val):
    try:
        return int(val)
    except ValueError, exc:
        return None

class search_error_wrapper(object):
    "Used as decorator to display  errors in the search calls"
    def __init__(self, f):
        self.f = f
        
    def __call__(self, *args, **kwds):
        try:
            # first parameter is the request
            value = self.f(*args, **kwds)
            return value
        except Exception,exc:
            request = kwds.get('request') or args[0]
            messages.error(request, "Search error: %s" % exc)
            return []
            
@search_error_wrapper         
def search_results(request, text, subset=None):
    "Returns search results for a text"
    text = text.strip()[:200]
    
    if not text:
        return []
    
    ix = index.open_dir(settings.WHOOSH_INDEX)
    searcher = ix.searcher()
    parser   = MultifieldParser(["content"], schema=ix.schema)
    #parser.remove_plugin_class(WildcardPlugin)
    query    = parser.parse(text)
    results  = searcher.search(query, limit=200)
    results.formatter.maxchars = 350
    results = map(decorate, results)
    if subset:
        results = filter(lambda r:r['type']in subset, results)
    searcher.close()
    ix.close()

    return results

def decorate(res):
    "Decorates search results with highlights"
    content = res.highlights('content')
    return html.Params(title=res['title'], pid=res['pid'], content=content, type=res['type'])

def get_subset(word):
    "Returns a set of post types based on a word"
    if word == 'all':
        subset = []
    elif word == 'top':
        subset = set( (POST_MAP[key] for key in POST_TOPLEVEL) )
    else:
        subset = [ word ]
    return subset

def main(request):
    "Main search"
    counts = request.session.get(SESSION_POST_COUNT, {})
    
    q = request.GET.get('q','') # query
    t = request.GET.get('t','all')  # type
    
    params = html.Params(tab='search', q=q, sort='')
    subset = get_subset(t)
   
    if params.q:
        form = SearchForm(request.GET)
        res  = search_results(request=request, text=params.q, subset=subset)
        objects = Post.objects.filter(id__in=[r['pid'] for r in res]).exclude(status=POST_DELETED)
        for object, r in zip(objects, res):
            object.context = r['content']
        size = len(res)
        #messages.info(request, 'Searched results for: %s found %d results' % (params.q, size))
    else:
        form = SearchForm()
        res  = []
        objects = []

    page = get_page(request, objects, per_page=10)
    return html.template(request, name='search.html', page=page, params=params, counts=counts, form=form)

# number of terms extracted during a more like this query
NUM_TERMS = 10
TOP_COUNT = 20

@search_error_wrapper  
def more_like_this(request, pid):
    ix = index.open_dir(settings.WHOOSH_INDEX)
    searcher = ix.searcher()
    qp = QueryParser("pid", schema=ix.schema)
    qq = qp.parse(pid)
    doc = searcher.search(qq)
     
    first = doc[0]
    title = "%s: %s" % (first['type'], first['title'])

    res = first.more_like_this("content", numterms=NUM_TERMS)
    res = map(decorate, res)
    ix.close()
    
    messages.info(request, 'Posts similar to <b>%s</b>' % title)

    return res

import time
def print_timing(func):
    def wrapper(*args, **kwds):
        start = time.time()
        res = func(*args, **kwds)
        end = time.time() 
        diff = (end - start ) * 1000
        print '**** function %s in %0.3f ms' % (func.func_name, diff)
        return res
    return wrapper

@print_timing
def more(request, pid):
    counts = request.session.get(SESSION_POST_COUNT, {})
    form = SearchForm()
    params = html.Params(tab='search', q="", sort='')
    res = more_like_this(request=request, pid=pid)
    page = get_page(request, res, per_page=10)
    return html.template(request, name='search.html', page=page, params=params, counts=counts, form=form)

def update(post, handler=None):
    "Adds/updates a post to the index"
    
    if handler:
        writer = handler
    else:
        ix = index.open_dir(settings.WHOOSH_INDEX)
        writer = ix.writer()
        
    pid = post.id
    content = unicode(post.content)
    title  = unicode(post.title)
    type = unicode(post.get_type_display())
    if post.type in POST_TOPLEVEL:
        content = unicode(post.title + " " + post.content)
    else:
        content = unicode(post.content)

    writer.update_document(pid=pid, content=content, title=title, type=type)
    
    # only commit if this was opened here
    if not handler:
        writer.commit()

def add_batch(posts):
    ix = index.open_dir(settings.WHOOSH_INDEX)
    wr = ix.writer(limitmb=50)
    for post in posts:
        update(post, handler=wr)
    wr.commit()
    ix.close()

def reset_index():
    "Resets the index"
    info("reset index at %s" % settings.WHOOSH_INDEX)
    ix = index.create_in(settings.WHOOSH_INDEX, SCHEMA)
    posts = get_posts(False)
    # reset the status of all posts
    posts.update(changed=True)

def get_posts(flag=True):
    from main.server import models
    return models.Post.objects.filter(changed=flag).exclude(type=POST_COMMENT)

def full_index():
    """
    Runs indexing on all changed posts
    """
   
    import glob

    path = "%s/*" % settings.WHOOSH_INDEX
    if not glob.glob(path):
       reset_index()

    # get a count of posts
    count = get_posts().count()
    info("indexing %s posts" % count)
   
    # index the posts
    posts = get_posts().all()
    add_batch(posts)
    
    # update the attribute
    posts = get_posts().update(changed=False)
        
if __name__ == '__main__':
    import doctest, optparse
   
    parser = optparse.OptionParser()
    parser.add_option("--reset", dest="reset",  action="store_true", default=False, help="resets the index")

    (opts, args) = parser.parse_args()
    if opts.reset:
        reset_index()

    full_index()