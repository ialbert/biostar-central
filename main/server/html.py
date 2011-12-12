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

import markdown, pygments
from pygments import highlight, lexers, formatters
from pygments.formatters import HtmlFormatter
from itertools import groupby

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
    """
    Generates html from a markdown text.
    
    Also custom renders codeblocks via pygments
    
    >>> generate("ABCD")
    """
    
    if not text:
        return ""

    md = markdown.Markdown(
        safe_mode=True,
    )
    md.html_replacement_text = "[?]"

    # split the text into lines
    lines = text.splitlines()
    
    # this is to save some effort on single line content
    if len(lines) == 1:
        html = md.convert(lines[0])
        return html

    # function to detect code starts
    func  = lambda line: line.startswith(' ' * 4) or line.startswith("\t")
    pairs = zip(map(func, lines), lines) 

    # mark single-newlines surrounded by code as code
    pairs2 = []
    for i, p in enumerate(pairs):
        flag, line = p
        if i != 0 and i < len(pairs) - 1  and \
          line == '' and pairs[i - 1][0] and pairs[i + 1][0]:
            pairs2.append( (True, '') )
        else:
            pairs2.append( (flag, line) )

    # grouping function, group by the codeblock condition
    func   = lambda pair: pair[0]
    blocks = groupby(pairs2, func)

    # tranform to continus text within each block
    groups = []
    for flag, group in blocks:
        block = '\n'.join( g[1] for g in group)
        groups.append( (flag, block) )

    # markup each block as needed
    out, code = [], {}
    for flag, block in groups:
        if flag:
            # this is codeblock
            try:
                # guess syntax highlighting
                lexer = lexers.guess_lexer(block)
            except pygments.util.ClassNotFound, exc:
                # unable to detect language fall back to Python
                lexer = lexers.PythonLexer()
            body = pygments.highlight(block, lexer, formatters.HtmlFormatter())
            uuid = hashlib.md5( str(random.getrandbits(256))).hexdigest()
            code[uuid] = body
            out.append(uuid)
        else:
            # regular markdown
            out.append( block )

    # markdown needs to run on the entire body to handle references
    text = '\n'.join(out)
    html = md.convert(text)
    
    # this is required to put the codes back into the markdown
    for uuid, value in code.items():
        html = html.replace(uuid, value)

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
        self.__dict__.update(kwds)
        
    def parse(self, request):
        if 'q' in request.GET:
            self.q = request['q']
            self.setr('Searching for %s' % self.q)
        else:
            self.q = None

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
