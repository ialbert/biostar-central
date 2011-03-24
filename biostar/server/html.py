"""
Html utility functions.
"""
import string, mimetypes, os
from django.template import Context, loader
from django.core.servers.basehttp import FileWrapper
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render_to_response
from django.core.context_processors import csrf


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

def redirect(url):
    "Redirects to a url"
    return HttpResponseRedirect(url)

def template(request, name, mimetype=None, **kwd):
    """Renders a template and returns it as an http response"""
    user = request.user        
    messages = user.get_and_delete_messages()
    
    # this collects the dictionary from which the context will be built
    cx = dict(messages=messages, user=user)
    cx.update(kwd)
    cx.update(csrf(request))
    
    context  = Context(cx)
    template = loader.get_template(name)
    page = template.render(context)

    return render_to_response(name, context)

def test():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    test()