__author__ = 'ialbert'
import json, traceback, logging
from braces.views import JSONResponseMixin
from biostar.apps.posts.models import Post
from biostar.apps.users.models import User
from django.views.generic import View
from django.shortcuts import render_to_response, render
from django.http import HttpResponse, HttpResponseRedirect, HttpResponsePermanentRedirect, Http404
from functools import partial

def json_response(adict, **kwd):
    """Returns a http response in JSON format from a dictionary"""
    return HttpResponse(json.dumps(adict), **kwd)


logger = logging.getLogger(__name__)

def ajax_msg(msg, status):
    return json_response(dict(status=status, msg=msg))

ajax_success = partial(ajax_msg, status='success')
ajax_error   = partial(ajax_msg, status='error')


class ajax_error_wrapper(object):
    """
    Used as decorator to trap/display  errors in the ajax calls
    """
    def __init__(self, f):
        self.f = f

    def __call__(self, request):
        try:
            if request.method != 'POST':
                return ajax_error('POST method must be used.')

            if not request.user.is_authenticated():
                return ajax_error('You must be logged in to do that')

            value = self.f(request)
            return value
        except Exception,exc:
            traceback.print_exc()
            return ajax_error('Error: %s' % exc)

@ajax_error_wrapper
def vote(request):
    "Handles all voting on posts"

    print(request.POST)
    return ajax_success("OK")