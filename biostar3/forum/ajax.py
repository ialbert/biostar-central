
import json, traceback, logging
from django.conf import settings
from django.shortcuts import render, redirect, render_to_response
from django.http import HttpResponse
from functools import partial
from django.db import transaction
from django.db.models import Q, F
from .models import Post, User, Vote, GroupPerm
from . import auth
from django.contrib import messages
from functools import partial

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

def add_comment(request, pk):
    template_name = "post_comment_add.html"
    context = dict(pk=pk)
    return render(request, template_name, context)

@auth.post_edit
def post_moderate(request, pk, post=None, user=None):
    """
    Authorization can get complicated due to its branching nature.

    To simply the paths we do not nest conditionals or
    write if/elif constructs.

    Instead match conditions and return from them as soon as possible.
    """
    template_name = "post_moderate.html"
    back = redirect(post.get_absolute_url())

    # These actions need to match the templates.
    CLOSE, DELETE, DUPLICATE, REPARENT, RESTORE = "close", "delete", "duplicate", "reparent", "restore"
    MOVE_TO_ANSWER, MOVE_TO_COMMENT, MOVE_TO_QUESTION = "move_to_answer", "move_to_comment",  "move_to_question"

    # Simplify a few functions
    error = partial(messages.error, request)
    info = partial(messages.info, request)
    select = Post.objects.filter(pk=post.id)

    if request.method != "POST":
        # Any other request gets a template.
        context = dict(post=post)
        return render(request, template_name, context)

    # Get the requested actions.
    action = request.POST.get("action")
    fields = request.POST.getlist("reason")
    reason = ''.join(fields)

    is_author = (post.author == user)

    # Get the parent id.
    parent_id = request.POST.get("parent_id", "0")
    parent_id = auth.safe_int(parent_id)
    parent = Post.objects.filter(pk=parent_id).first()

    if not auth.can_moderate_post(user=user, post=post):
        error("You have insufficient permissions to moderate that post")
        return back

    if parent_id and not parent:
        error("Invalid parent id.")
        return back

    if parent == post:
        error("Re-parenting to the same post.")
        return back

    if parent and parent.root != post.root:
        error("May only reparent in the same thread.")
        return back

    if action in (CLOSE, DELETE, DUPLICATE) and not reason:
        error("Must specify a reason for a close/delete.")
        return back

    if action == RESTORE and not user.is_moderator:
        error("Only moderators may restore posts")
        return back

    if not post.is_toplevel and action in (CLOSE, DUPLICATE):
        error("Only top level posts may be closed")
        return back

    if post.is_toplevel and action in (MOVE_TO_ANSWER, MOVE_TO_COMMENT, MOVE_TO_QUESTION):
        error("Top level post may not be moved around.")
        return back

    # At this point the parameters are correct.

    # Post creation shortcut
    create = partial(auth.create_content_post,
                     content=reason, parent=post, user=user)

    if action == CLOSE:
        info("Post closed")
        select.update(status=Post.CLOSED)
        create(post_type=Post.COMMENT)
        return back

    if action == RESTORE:
        info("Post restored")
        select.update(status=Post.OPEN)
        return back

    if action == DUPLICATE:
        info("Post closed as duplicate")
        select.update(status=Post.CLOSED)
        create(post_type=Post.ANSWER)
        return back

    if action == DELETE:
        info("Post deleted")
        select.update(status=Post.DELETED)
        create(post_type=Post.COMMENT)
        if is_author and (post.reply_count == 0):
            post.delete()
        return back

    if action == MOVE_TO_ANSWER:
        info("Post moved to an answer")
        select.update(post_type=Post.ANSWER)
        post.root.set_reply_count()
        return back

    if action == MOVE_TO_COMMENT:
        info("Post moved to a comment")
        select.update(type=Post.COMMENT)
        post.root.set_reply_count()
        return back

    if action == MOVE_TO_QUESTION:
        info("Post moved to a new question")
        pass

    if action == REPARENT:
        info("Post reparented")
        select.update(parent=parent)
        return back

    return back

@auth.valid_user
def user_moderate(request, pk, target=None):
    template_name = "user_moderate.html"

    back = redirect(target.get_absolute_url())

    # Simplify a few functions
    error = partial(messages.error, request)
    info = partial(messages.info, request)
    select = User.objects.filter(pk=target.id)

    if request.method != "POST":
        context = dict(target=target)
        return render(request, template_name, context)

    info("Moderation completed")

    return back