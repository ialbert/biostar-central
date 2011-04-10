"""
Html utility functions.
"""
import re, string, mimetypes, os, json
from django.template import RequestContext, loader
from django.core.servers.basehttp import FileWrapper
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render_to_response
from django.core.context_processors import csrf
from BeautifulSoup import BeautifulSoup, Comment

import markdown, pygments
from pygments import highlight, lexers, formatters
from pygments.formatters import HtmlFormatter
from itertools import groupby

def generate(text):
    """
    Generates html from a markdown text
    
    >>> generate("ABCD")
    """

    # split the text into lines
    lines = text.splitlines()
    
    # function to detect code starts
    func  = lambda line: line.startswith(' ' * 4) or line.startswith("\t")
    pairs = zip(map(func, lines), lines) 

    # grouping function, group by the codeblock condition
    func   = lambda pair: pair[0]
    blocks = groupby(pairs, func)

    # tranform to continus text within each block
    groups = []
    for flag, group in blocks:
        block = '\n'.join( g[1] for g in group)
        groups.append( (flag, block) )

    # markup each block as needed
    out = []
    for flag, block in groups:
        if flag:
            # this is codeblock
            try:
                # guess syntax highlighting
                lexer = lexers.guess_lexer(block)
            except pygments.util.ClassNotFound, exc:
                # unable to detect language fall back to Python
                lexer = lexers.PythonLexer()
            body  = pygments.highlight(block, lexer, formatters.HtmlFormatter())
        else:
            # regular markdown
            body = markdown.markdown(block, safe_mode=True)
        out.append( body )

    html = '\n'.join(out)
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
    Represents incoming parameters. Keyword arguments
    will be defaults.

    >>> p = Params(a=1, b=2, c=3, incoming=dict(c=4))
    >>> p.a, p.b, p.c
    (1, 2, 4)
    """
    def __init__(self, incoming={}, **kwds):
        self.__dict__.update(kwds)
        self.__dict__.update(incoming)

    def update(self, data):
        self.__dict__.update(data)

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self):
        return 'Params: %s' % self.__dict__

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

def test():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    test()