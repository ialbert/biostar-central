
import json, traceback, logging
from django.conf import settings
from django.http import HttpResponse
from functools import partial
from django.db import transaction
from django.db.models import Q, F
from .models import Post, User, Vote

logger = logging.getLogger("biostar")

def json_response(data, **kwd):
    """
    Returns data as a http response in JSON format
    """
    return HttpResponse(json.dumps(data), **kwd)


def ajax_msg(msg, status, **kwargs):
    """
    Returns a json response that encodes a dictionary with
    at least two keys: 'msg' and 'status' plus optional additional data.

    All Biostar responses will follow this format. The status will indicate 'success' or 'error'.
    The message is to describe the status. The rest of the content is data.
    """
    data = dict(status=status, msg=msg)
    data.update(kwargs)
    return json_response(data)

# Shortcuts for success and error handlers.
ajax_error = partial(ajax_msg, status='error')
ajax_success = partial(ajax_msg, status='success')

class ajax_error_wrapper(object):
    """
    Used as decorator to trap errors in the ajax calls
    """

    def __init__(self, f):
        self.f = f

    def __call__(self, request):
        try:
            if request.method != 'POST':
                return ajax_error('POST method must be used.')

            if not request.user.is_authenticated():
                return ajax_error('You must be logged in to access the site')

            value = self.f(request)
            return value
        except Exception, exc:
            if settings.DEBUG:
                traceback.print_exc()
            logger.error(exc)
            return ajax_error('Error: %s' % exc)

@ajax_error_wrapper
def vote_handler(request):
    """
    Handles voting on posts
    """
    msg = "Success!"
    return ajax_success(msg)
