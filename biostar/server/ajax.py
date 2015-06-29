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


def ajax_msg(msg, status, **kwargs):
    payload = dict(status=status, msg=msg)
    payload.update(kwargs)
    return json_response(payload)


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


POST_TYPE_MAP = dict(vote=Vote.UP, bookmark=Vote.BOOKMARK, accept=Vote.ACCEPT)

@transaction.atomic
def perform_vote(post, user, vote_type):

    # Only maintain one vote for each user/post pair.
    votes = Vote.objects.filter(author=user, post=post, type=vote_type)
    if votes:
        vote = votes[0]
        msg = "%s removed" % vote.get_type_display()
        change = -1
    else:
        change = +1
        vote = Vote.objects.create(author=user, post=post, type=vote_type)
        msg = "%s added" % vote.get_type_display()

    if post.author != user:
        # Update the user reputation only if the author is different.
        User.objects.filter(pk=post.author.id).update(score=F('score') + change)

    # The thread score represents all votes in a thread
    Post.objects.filter(pk=post.root_id).update(thread_score=F('thread_score') + change)

    if vote.type == Vote.BOOKMARK:
        # Apply the vote
        Post.objects.filter(pk=post.id).update(book_count=F('book_count') + change, vote_count=F('vote_count') + change)
        Post.objects.filter(pk=post.id).update(subs_count=F('subs_count') + change)
        Post.objects.filter(pk=post.root_id).update(subs_count=F('subs_count') + change)

    elif vote_type == Vote.ACCEPT:
        if change > 0:
            # There does not seem to be a negation operator for F objects.
            Post.objects.filter(pk=post.id).update(vote_count=F('vote_count') + change, has_accepted=True)
            Post.objects.filter(pk=post.root_id).update(has_accepted=True)
        else:
            Post.objects.filter(pk=post.id).update(vote_count=F('vote_count') + change, has_accepted=False)

            # Only set root as not accepted if there are no accepted siblings
            if Post.objects.exclude(pk=post.root_id).filter(root_id=post.root_id, has_accepted=True).count() == 0:
                Post.objects.filter(pk=post.root_id).update(has_accepted=False)
    else:
        Post.objects.filter(pk=post.id).update(vote_count=F('vote_count') + change)

    # Clear old votes.
    if votes:
        votes.delete()

    return msg




@ajax_error_wrapper
def vote_handler(request):
    "Handles all voting on posts"


    user = request.user
    vote_type = request.POST['vote_type']
    vote_type = POST_TYPE_MAP[vote_type]
    post_id = request.POST['post_id']

    # Check the post that is voted on.
    post = Post.objects.get(pk=post_id)

    if post.author == user and vote_type == Vote.UP:
        return ajax_error("You can't upvote your own post.")

    #if post.author == user and vote_type == Vote.ACCEPT:
    #    return ajax_error("You can't accept your own post.")

    if post.root.author != user and vote_type == Vote.ACCEPT:
        return ajax_error("Only the person asking the question may accept this answer.")

    with transaction.atomic():
        msg = perform_vote(post=post, user=user, vote_type=vote_type)

    return ajax_success(msg)