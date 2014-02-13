__author__ = 'ialbert'
import json, traceback, logging
from braces.views import JSONResponseMixin
from biostar.apps.posts.models import Post, Vote
from biostar.apps.users.models import User
from django.views.generic import View
from django.shortcuts import render_to_response, render
from django.http import HttpResponse, HttpResponseRedirect, HttpResponsePermanentRedirect, Http404
from functools import partial
from django.db import transaction
from django.db.models import Q, F


def json_response(adict, **kwd):
    """Returns a http response in JSON format from a dictionary"""
    return HttpResponse(json.dumps(adict), **kwd)


logger = logging.getLogger(__name__)


def ajax_msg(msg, status):
    return json_response(dict(status=status, msg=msg))


ajax_success = partial(ajax_msg, status='success')
ajax_error = partial(ajax_msg, status='error')


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
        except Exception, exc:
            traceback.print_exc()
            return ajax_error('Error: %s' % exc)


@ajax_error_wrapper
def vote(request):
    "Handles all voting on posts"

    user = request.user
    post_type = request.POST['post_type']
    post_id = request.POST['post_id']

    # Check the post that is voted on.
    post = Post.objects.get(pk=post_id)

    if post.author == user and post_type == "vote":
        return ajax_error("May not vote on your own posts")

    with transaction.atomic():
        # Only maintain one vote for each user/post pair.
        votes = Vote.objects.filter(author=user, post=post)
        if votes:
            change = -1
            votes.delete()
        else:
            change = +1
            Vote.objects.create(author=user, post=post, type=Vote.UP)

        # Update the user score.
        User.objects.filter(pk=post.author.id).update(score=F('score') + change)

    return ajax_success("OK")