"""
Html utility functions.
"""
import re, string, mimetypes, os, json, random, hashlib,  unittest
from django.template import RequestContext, loader
from django.core.servers.basehttp import FileWrapper
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render_to_response
from django.core.context_processors import csrf
from django.core.paginator import Paginator, InvalidPage, EmptyPage

from BeautifulSoup import BeautifulSoup, Comment

import markdown
from docutils import core

from itertools import groupby

# safe string transformation
import string
VALID = set(string.ascii_letters + string.digits + "-_ ")

def ascii(text):
    global VALID
    return filter(lambda x: x in VALID, text)

def get_page(request, obj_list, per_page=25):
    "A generic paginator"

    paginator = Paginator(obj_list, per_page) 
    try:
        pid = int(request.GET.get('page', '1'))
    except ValueError:
        pid = 1

    try:
        page = paginator.page(pid)
    except (EmptyPage, InvalidPage):
        page = paginator.page(paginator.num_pages)
    
    return page

def nuke(text):
    """
    This function is not the main sanitizer,
    is used mainly as an extra precaution to preemtively
    delete markup from markdown content.
    """
    text = text.replace("<","&lt;")
    text = text.replace(">","&gt;")
    text = text.replace("\"","&quot;")
    text = text.replace("&","&amp;")
    return text

def generate(text):
    if not text:
        return ""
    text = text.rstrip()
    if text.startswith('##rest'):
        text = text[6:].strip()
        rest = core.publish_parts(text ,writer_name='html')
        html = rest.get('html_body','[rest error]')
    else:
        md = markdown.Markdown( safe_mode=True )
        md.html_replacement_text = "[?]"
        html = md.convert(text)
    return html

ALLOWED_TAGS = "strong span:class br ol ul li a:href img:src pre code blockquote p em"
def sanitize(value, allowed_tags=ALLOWED_TAGS):
    """
    From http://djangosnippets.org/snippets/1655/

    Argument should be in form 'tag2:attr1:attr2 tag2:attr1 tag3', where tags
    are allowed HTML tags, and attrs are the allowed attributes for that tag.
    """
    js_regex = re.compile(r'[\s]*(&#x.{1,7})?'.join(list('javascript')))
    allowed_tags = [tag.split(':') for tag in allowed_tags.split()]
    allowed_tags = dict((tag[0], tag[1:]) for tag in allowed_tags)

    soup = BeautifulSoup(value)
    for comment in soup.findAll(text=lambda text: isinstance(text, Comment)):
        comment.extract()

    for tag in soup.findAll(True):
        if tag.name not in allowed_tags:
            tag.hidden = True
        else:
            tag.attrs = [(attr, js_regex.sub('', val)) for attr, val in tag.attrs
                         if attr in allowed_tags[tag.name]]

    return soup.renderContents().decode('utf8')

class Params(object):
    """
    Represents incoming parameters. 
    Parameters with special meaning: q - search query, m - matching
    Keyword arguments
    will be defaults.

    >>> p = Params(a=1, b=2, c=3, incoming=dict(c=4))
    >>> p.a, p.b, p.c
    (1, 2, 4)
    """
    def __init__(self, **kwds):
        self.q = ''
        self.__dict__.update(kwds)
        
    def parse(self, request):
        self.q = request.GET.get('q', '')
        if self.q:
            self.setr('Searching for %s' % self.q)
    
    def get(self, key, default=''):
        return self.__dict__.get(key, default)
        
    def update(self, data):
        self.__dict__.update(data)

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self):
        return 'Params: %s' % self.__dict__
    
    def setr(self, text):
        setattr(self, 'remind', text)

    def getr(self,text):
        return getattr(self, 'remind')

def response(data, **kwd):
    """Returns a http response"""
    return HttpResponse(data, **kwd)
    
def json_response(adict, **kwd):
    """Returns a http response in JSON format from a dictionary"""
    return HttpResponse(json.dumps(adict), **kwd)

def redirect(url):
    "Redirects to a url"
    return HttpResponseRedirect(url)

def template(request, name, mimetype=None, **kwd):
    """Renders a template and returns it as an http response"""
    
    # parameters that will always be available for the template
    kwd['request'] = request
    return render_to_response(name, kwd, context_instance=RequestContext(request))
    
class HtmlTest(unittest.TestCase):
    def test_sanitize(self):
        "Testing HTML sanitization"
        text = sanitize('<a href="javascrip:something">A</a>', allowed_tags="b")
        self.assertEqual( text, u'A' )

        markup = generate("ABCD\n    CDE\n*A*")
        expect = '<p>ABCD\n<div class="highlight"><pre>    <span class="n">CDE</span>\n</pre></div>\n\n<em>A</em></p>'
        self.assertEqual( markup, expect )

        p = Params(a=1, b=2, c=3)
        self.assertEqual( (p.a, p.b, p.c), (1, 2, 3))

def suite():
    s = unittest.TestLoader().loadTestsFromTestCase(HtmlTest)
    return s
