
import json, traceback, logging
from django.conf import settings
from django.shortcuts import render, redirect, render_to_response
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
                return ajax_error('You must be logged to do that!')

            value = self.f(request)
            return value
        except Exception, exc:
            if settings.DEBUG:
                traceback.print_exc()
            logger.error(exc)
            return ajax_error('Error: %s' % exc)


VOTE_TYPE_MAP = dict(upvote=Vote.UP, bookmark=Vote.BOOKMARK, accept=Vote.ACCEPT)

@transaction.atomic
def perform_vote(post, user, vote_type):
    # Only maintain one vote for each user/post pair.
    # Normally one should never end up with more than one vote per user/post but
    # some databases don't enforce this and/or concurrency problems may lead to this.
    votes = Vote.objects.filter(author=user, post=post, type=vote_type)

    if votes:
        vote = votes[0]
        msg = "%s removed" % vote.get_type_display()
        change = -len(votes)
    else:
        change = +1
        vote = Vote.objects.create(author=user, post=post, type=vote_type)
        msg = "%s added" % vote.get_type_display()

    if post.author != user:
        # Update the user reputation only if the author is different.
        User.objects.filter(pk=post.author.id).update(score=F('score') + change)

    # The thread score represents all votes in a thread.
    Post.objects.filter(pk=post.root_id).update(thread_score=F('thread_score') + change)

    if vote.type == Vote.BOOKMARK:
        # Bookmarks also alter attributes on the root post.
        Post.objects.filter(pk=post.id).update(book_count=F('book_count') + change, vote_count=F('vote_count') + change)
        Post.objects.filter(pk=post.id).update(subs_count=F('subs_count') + change)
        Post.objects.filter(pk=post.root_id).update(subs_count=F('subs_count') + change)

    elif vote_type == Vote.ACCEPT:
        # There does not seem to be a negation operator for F objects hence
        # this ends up a little more complicated than it should be.
        if change > 0:

            Post.objects.filter(pk=post.id).update(vote_count=F('vote_count') + change, has_accepted=True)
            Post.objects.filter(pk=post.root_id).update(has_accepted=True)
        else:
            Post.objects.filter(pk=post.id).update(vote_count=F('vote_count') + change, has_accepted=False)
            Post.objects.filter(pk=post.root_id).update(has_accepted=False)
    else:
        Post.objects.filter(pk=post.id).update(vote_count=F('vote_count') + change)

    # Clear old votes.
    if votes:
        votes.delete()

    return msg


@ajax_error_wrapper
def vote_handler(request):
    """
    Handles voting on posts
    """
    global VOTE_TYPE_MAP

    user = request.user
    vote_type = request.POST['vote_type']
    vote_type = VOTE_TYPE_MAP[vote_type]
    post_id = request.POST['post_id']

    # Check the post that is voted on.
    post = Post.objects.get(pk=post_id)

    if post.author == user and vote_type == Vote.UP:
        return ajax_error("You can't upvote your own post.")

    msg = perform_vote(post=post, user=user, vote_type=vote_type)

    return ajax_success(msg)

TEMPLATE_MAPPER = dict(
    comment_panel='widgets/comment_panel.html',
)
def load_html(request, name, pk):
    global TEMPLATE_MAPPER
    default = TEMPLATE_MAPPER['comment_panel']
    template_name = TEMPLATE_MAPPER.get(name, default)
    context = dict(pk=pk)
    return render(request, template_name, context)
