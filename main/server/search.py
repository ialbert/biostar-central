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
from whoosh.analysis import StemmingAnalyzer
from django.contrib import messages
from itertools import *

stem = StemmingAnalyzer()
SCHEMA = fields.Schema(
    title   = fields.TEXT(analyzer=stem, stored=True),
    content = fields.TEXT(analyzer=stem, stored=True), 
    type    = fields.TEXT(stored=True),
    pid     = fields.NUMERIC(stored=True, unique=True),
    uid     = fields.NUMERIC(stored=True),
    author  = fields.TEXT(stored=True) 
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
    q = forms.CharField(max_length=30,  initial="", widget=forms.TextInput(attrs={'size':'50'}))   
    t = forms.ChoiceField(choices=choices, required=False)

def safe_int(val):
    try:
        return int(val)
    except ValueError, exc:
        return None
    
def search_query(request, text, subset=None):
    text = text.strip()[:200]
    
    if not text:
        return []
    
    try:
        ix = index.open_dir(settings.WHOOSH_INDEX)
        searcher = ix.searcher()
        parser   = MultifieldParser(["title", "content"], schema=ix.schema)
        #parser.remove_plugin_class(WildcardPlugin)
        query    = parser.parse(text)
        results  = searcher.search(query, limit=200)
        results.formatter.maxchars = 350
        results = map(decorate, results)
        if subset:
            results = filter(lambda r:r['type']in subset, results)
        searcher.close()
        ix.close()
    except Exception, exc:
        messages.error(request, "Search error: %s" % exc)
        results = []
   
    return results

def decorate(res):
    content = res.highlights('content')
    return html.Params(title=res['title'], uid=res['uid'],
                       pid=res['pid'], content=content, type=res['type'])

def main(request):
    
    counts = request.session.get(SESSION_POST_COUNT, {})
    
    q = request.GET.get('q','') # query
    t = request.GET.get('t','all')  # type
    
    params = html.Params(tab='search', q=q)

    if t == 'all':
        subset = None
    elif t == 'top':
        subset = set( (POST_MAP[key] for key in POST_TOPLEVEL) )
    else:
        subset = [ t ]

    if params.q:
        form = SearchForm(request.GET)
        res  = search_query(request=request, text=params.q, subset=subset)
        size = len(res)
        messages.info(request, 'Searched results for: %s found %d results' % (params.q, size))
    else:
        form = SearchForm()
        res  = []
    
    page = get_page(request, res, per_page=10)
    return html.template(request, name='search.html', page=page, params=params, counts=counts, form=form)

def more(request, pid):
    counts = request.session.get(SESSION_POST_COUNT, {})
    form = SearchForm()
    
    params = html.Params(tab='search', q=pid)

    ix = index.open_dir(settings.WHOOSH_INDEX)
    searcher = ix.searcher()
    qp = QueryParser("pid", schema=ix.schema)
    qq = qp.parse(pid)
    rr = searcher.search(qq)
    
    if not rr:
        messages.error(request, 'Server settings problem - this post is not indexed! Please report.')
        res = []
    else:
        first = rr[0]
        messages.info(request, 'Searching for posts similar to: %s' % first['title'])
    
        subset = set( (POST_MAP[key] for key in POST_TOPLEVEL) )
        
        res = first.more_like_this("content")
        res = filter(lambda x: x['type'] in subset, res)
        res = map(decorate, res)
        
    page = get_page(request, res, per_page=10)
    return html.template(request, name='search.html', page=page, params=params, counts=counts, form=form)

def update(post, created, handler=None):
    "Adds/updates a post to the index"
    
    if handler:
        writer = handler
    else:
        ix = index.open_dir(settings.WHOOSH_INDEX)
        writer = ix.writer()
        
    # set the author name
    content = post.content
    title   = post.title
    author  = unicode(post.author.profile.display_name)
    uid     = post.author.id
    content = unicode(content)
    title   = unicode(title)
    type    = unicode(post.get_type_display())
    
    if created:                     
        writer.add_document(content=content, pid=post.id, author=author, title=title, uid=uid, type=type)
    else:
        writer.update_document(content=content, pid=post.id, author=author, title=title, uid=uid, type=type)
    
    # only commit if this was opened here
    if not handler:
        writer.commit()

def add_batch(posts):
    gc.collect()
    ix = index.open_dir(settings.WHOOSH_INDEX)
    wr = ix.writer(limitmb=10)
    for post in posts:
        update(post, created=True, handler=wr)
    wr.commit()
    ix.close()
    del wr, ix
    gc.collect()
    
def full_index():
    "Runs a full indexing on all posts"
    from main.server import models
    
    ix = index.create_in(settings.WHOOSH_INDEX, SCHEMA)
    count = models.Post.objects.all().count()
    info("found %s posts" % count)
    posts = []
    STEP  = 1000
    for lo in xrange(0, count, STEP):
        hi = min( (count, lo + STEP) )
        info("whoosh indexing %s posts" % hi)
        
        posts = list(models.Post.objects.all()[lo:hi])
        add_batch(posts)
        del posts
        gc.collect()
        
if __name__ == '__main__':
    full_index()