from __future__ import print_function, unicode_literals, absolute_import, division
import json, traceback, logging, pyzmail, ftfy
from django.conf import settings
from django.shortcuts import render, redirect, render_to_response
from django.http import HttpResponse
from functools import partial
from django.db import transaction
from django.db.models import Q, F
from .models import Post, User, Vote, ReplyToken
from . import auth
from django.contrib import messages
from functools import partial
from django.views.decorators.csrf import csrf_exempt
from biostar3.utils import email_reply_parser
from biostar3.utils.compat import *
from django.http import HttpResponse
from django.utils.encoding import smart_text

logger = logging.getLogger("biostar")


def json_response(data, **kwd):
    """
    Returns data as a http response in JSON format
    """
    return HttpResponse(json.dumps(data), **kwd)


NEW_POST, REPLY_POST = "new", "reply"

@csrf_exempt
def email_handler(request):
    """
    Handles email based interactions. The emails need to be posted
    in mbox type format.
    """

    # The feature may be turned off.
    if not settings.ALLOW_EMAIL_REPLY:
        resp = dict(status="error", msg="email reply not allowed")
        return json_response(resp)

    # The key to validate the post.
    key = request.POST.get("key")

    if key != settings.EMAIL_HANDLER_SECRET_KEY:
        resp = dict(status="error", msg="key does not match")
        return json_response(resp)

    # This will contain the incoming email message.
    content = request.POST.get("content")

    try:
        msg = pyzmail.PyzMessage.factory(content)
    except Exception as exc:
        logger.error(exc)
        content = content.encode('utf8', errors='ignore')
        msg = pyzmail.PyzMessage.factory(content)

    # Extract the address from the address tuples.
    # The expected format is: action+token+group@site.com
    try:
        subject = msg.get_subject()
        address = msg.get_addresses('to')[0][1]
        target, site = address.split('@')
    except Exception as exc:
        resp = dict(status="error", msg="Invalid address={}".format(exc))
        return json_response(resp)

    # Parse the email address.
    parts = target.split('+')
    if len(parts) != 3:
        resp = dict(status="error", msg="Target address={} has incorrect format.".format(address))
        return json_response(resp)

    # Assign the various parts of the address.
    action, token, pattern = parts

    # Find the email content and its encoding
    part = msg.text_part or msg.html_part
    enc = part.charset
    text = part.get_payload().decode(enc)

    # Should we remove quote text.
    if settings.EMAIL_REPLY_REMOVE_QUOTED_TEXT:
        text = email_reply_parser.EmailReplyParser.parse_reply(text)

    # A new post is to be generated.
    if action == NEW_POST:

        # Which group to post to.
        usergroup = UserGroup.objects.filter(domain=pattern).first()
        if not usergroup:
            resp = dict(status="error", msg="Group={} does not exist".format(pattern))
            return json_response(resp)

        # The token must match a users uuid.
        author = User.objects.filter(profile__uuid=token).first()
        if not author:
            resp = dict(status="error", msg="Author uuid={} does not exist".format(token))
            return json_response(resp)

        # Create the post in the usergroup
        data = dict(title=subject, content=text, tags="via email")
        auth.create_toplevel_post(data=data, user=author, group=usergroup)
        resp = dict(status="ok", msg="question created".format(token))
        return json_response(resp)

    # Post a reply to a reply token.
    if action == REPLY_POST:
        replytok = ReplyToken.objects.filter(token=token).first()
        if not replytok:
            resp = dict(status="ok", msg="invalid reply token={}".format(token))
            return json_response(resp)
        auth.create_content_post(content=text, parent=replytok.post, user=replytok.user)
        resp = dict(status="error", msg="post created")
        return json_response(resp)

    resp = dict(status="error", msg="action not recognized")
    return json_response(resp)

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
        except Exception as exc:
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
        change = -1
    else:
        change = +1
        vote = Vote.objects.create(author=user, post=post, type=vote_type)
        msg = "%s added" % vote.get_type_display()

    if post.author != user:
        # Update the user reputation only if the author is different.
        User.objects.filter(pk=post.author.id).update(score=F('score') + change * 10)

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

    if vote_type == Vote.ACCEPT:
        if post.type != Post.ANSWER:
            return ajax_error("Only answers may be accepted.")

        if post.root.author != user:
            return ajax_error("Only the author of the top post may accept the answer.")


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
    MOVE_TO_ANSWER, MOVE_TO_COMMENT, MOVE_TO_QUESTION = "move_to_answer", "move_to_comment", "move_to_question"

    moderation_actions = post.moderation_actions(user)

    # Simplify a few functions
    error = partial(messages.error, request)
    info = partial(messages.info, request)
    select = Post.objects.filter(pk=post.id)

    if request.method != "POST":
        # Any other request gets a template.
        context = dict(post=post, moderation_actions=moderation_actions)
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

    if not auth.can_moderate_post(request=request, user=user, post=post):
        error("You have insufficient permissions to moderate that post")
        return back

    if not action in moderation_actions:
        error("Invalid moderation action for this post")
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
    user = request.user
    template_name = "user_moderate.html"

    back = redirect(target.get_absolute_url())

    # Simplify a few functions
    error = partial(messages.error, request)
    info = partial(messages.info, request)
    select = User.objects.filter(pk=target.id)

    if request.method != "POST":
        context = dict(target=target)
        return render(request, template_name, context)

    if not auth.can_moderate_user(request=request, user=request.user, target=target):
        error("You may not moderate that user")
        return back

    if user == target:
        error("You may not moderate yourself")
        return back

    action = request.POST.get("action")

    if not action:
        error("You must select a moderation action")
        return back

    SUSPEND, BAN, REINSTATE, MERGE = "suspend", "ban", "reinstate", "merge"

    if action == SUSPEND:
        select.update(status=User.SUSPENDED)
        info("User suspended")
        return back

    # Only staff may ban a user.
    if action == BAN and user.is_staff:
        select.update(status=User.BANNED, html="", content="")
        Post.objects.filter(author=target).delete()
        info("User banned")
        return back

    if action == REINSTATE:
        select.update(status=User.NEW_USER)
        info("User reinstated")
        return back

    info("Moderation completed")

    return back