"""
Indexes all post content
"""    
    
import shutil, os
from django import forms
from django.conf import settings
from main.server.const import *
from itertools import *
from main.server import html, const, formdef
from main.server.html import get_page
from whoosh import store, fields, index, highlight
from whoosh.qparser import QueryParser
from django.contrib import messages

# activate logging
import logging
logger = logging.getLogger(__name__)

SCHEMA = fields.Schema(
    content=fields.TEXT(stored=True), type=fields.TEXT(stored=True),
    pid=fields.NUMERIC(stored=True), uid=fields.NUMERIC(stored=True),
    title=fields.TEXT(stored=True), author=fields.TEXT(stored=True) )

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
    
class SearchForm(forms.Form):
    "A form representing a new question"
    q = forms.CharField(max_length=30,  initial="", widget=forms.TextInput(attrs={'size':'50'}))   
    type = forms.ChoiceField(choices= [ ("all", "All types") ] + list(POST_TYPES) , required=False)

def search_query(text):
    text = text.strip()[:200]
    if not text:
        return []
    ix = index.open_dir(settings.WHOOSH_INDEX)
    searcher = ix.searcher()
    parser   = QueryParser("content", ix.schema)
    query    = parser.parse(text)
    results  = searcher.search(query, limit=200)
    results.formatter.maxchars = 350
    return results

def decorate(res):
    
    content = res.highlights('content')
    return html.Params(title=res['title'], uid=res['uid'],
                       pid=res['pid'], content=content, type=res['type'])
       
def main(request):
    
    counts = request.session.get(SESSION_POST_COUNT, {})
    
    q = request.GET.get('q','')
    params = html.Params(tab='search', q=q)

    if params.q:
        form = SearchForm(request.GET)
        res  = search_query(params.q)
        res  = map(decorate, res)
        size = len(res)
        messages.info(request, 'Searched results for: %s found %d results' % (params.q, size))
        
    else:
        form = SearchForm()
        res  = []
    
    page = get_page(request, res, per_page=5)
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
        # updating
        writer.update_document(content=content, pid=post.id, author=author, title=title, uid=uid, type=type)
    
    # only commit if this was opened here
    if not handler:
        writer.commit()

def full_index():
    "Runs a full indexing on all posts"
    from main.server import models
    
    ix = index.create_in(settings.WHOOSH_INDEX, SCHEMA)
    wr = ix.writer()

    print "*** whoosh indexing %s posts" % models.Post.objects.all().count()
    for step, post in izip(count(1), models.Post.objects.all()):
        update(post, created=True, handler=wr)
        if step % 1000 == 0:
            wr.commit()
            wr = ix.writer()
    # final commit
    wr.commit()

if __name__ == '__main__':
    full_index()