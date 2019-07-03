from functools import wraps, partial
import logging
from ratelimit.decorators import ratelimit

from django.http import JsonResponse
from django.utils.decorators import available_attrs

from .const import *
from . import auth, util, forms, tasks, search
from .models import Post, Vote, Subscription


def ajax_msg(msg, status, **kwargs):
    payload = dict(status=status, msg=msg)
    payload.update(kwargs)
    return JsonResponse(payload)


logger = logging.getLevelName("biostar")
ajax_success = partial(ajax_msg, status='success')
ajax_error = partial(ajax_msg, status='error')


class ajax_error_wrapper:
    """
    Used as decorator to trap/display  errors in the ajax calls
    """

    def __init__(self, method):
        self.method = method

    def __call__(self, func, *args, **kwargs):

        @wraps(func, assigned=available_attrs(func))
        def _ajax_view(request, *args, **kwargs):
            if request.method != self.method:
                return ajax_error(f'{self.method} method must be used.')
            if not request.user.is_authenticated:
                return ajax_error('You must be logged in.')

            return func(request, *args, **kwargs)

        return _ajax_view


def ajax_test(request):
    """
    Creates a commment on a top level post.
    """
    msg="OK"
    print (f"HeRe= {request.POST} ")
    return ajax_error(msg=msg)


@ratelimit(key='ip', rate='500/h')
@ratelimit(key='ip', rate='25/m')
@ajax_error_wrapper(method="POST")
def ajax_vote(request):
    was_limited = getattr(request, 'limited', False)

    if was_limited:
        return ajax_error(msg="Too many votes from same IP address. Temporary ban.")

    user = request.user
    type_map = dict(upvote=Vote.UP, bookmark=Vote.BOOKMARK, accept=Vote.ACCEPT)

    vote_type = request.POST.get('vote_type')
    vote_type = type_map.get(vote_type)
    post_uid = request.POST.get('post_uid')

    # Check the post that is voted on.
    post = Post.objects.filter(uid=post_uid).first()

    if post.author == user and vote_type == Vote.UP:
        return ajax_error("You can not upvote your own post.")

    if post.author == user and vote_type == Vote.ACCEPT:
        return ajax_error("You can not accept your own post.")

    if post.root.author != user and vote_type == Vote.ACCEPT:
        return ajax_error("Only the person asking the question may accept this answer.")

    msg, vote, change = auth.apply_vote(post=post, user=user, vote_type=vote_type)

    return ajax_success(msg=msg, change=change)


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def ajax_subs(request):
    was_limited = getattr(request, 'limited', False)

    if was_limited:
        return ajax_error(msg="Too many votes from same IP address. Temporary ban.")

    type_map = dict(messages=Subscription.LOCAL_MESSAGE, email=Subscription.EMAIL_MESSAGE,
                    unfollow=Subscription.NO_MESSAGES)

    # Get the root and sub type.
    root_uid = request.POST.get('root_uid')
    sub_type = request.POST.get("sub_type")
    sub_type = type_map.get(sub_type, Subscription.NO_MESSAGES)
    user = request.user

    # Get the post that is subscribed to.
    root = Post.objects.filter(uid=root_uid).first()

    auth.create_subscription(post=root, user=user, sub_type=sub_type)

    return ajax_success(msg="Changed subscription.")


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def ajax_edit(request):
    """
    Edit post content using ajax.
    """
    uid = request.POST.get("post_uid")
    post = Post.objects.filter(uid=uid).first()

    if not post:
        return ajax_error(msg="Post does not exist")
    content = request.POST.get("content", post.content)
    length = len(content.replace(" ", ''))

    if length < forms.MIN_CONTENT:
        return ajax_error(msg=f"Too short, please add more than add more {forms.MIN_CONTENT} characters.")
    if length > forms.MAX_CONTENT:
        return ajax_error(msg=f"Too long, please add less than {forms.MAX_CONTENT} characters.")

    post.lastedit_user = request.user
    post.content = content
    post.save()

    # Note: returns html instead of JSON on success.
    # Used to switch content inplace.
    return ajax_success(msg=post.html)


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def ajax_search(request):

    query = request.POST.get('query', '')

    fields = ['content', 'tags', 'title']
    if query:
        results = search.search_index(query=query, fields=fields)

        print(results)
        1/0

    return



def validate_moderation(request, post):

    action = request.POST.get("action")
    dupe = request.POST.get("dupe")
    dupes = dupe.split(",")[:5]
    dupe_comment = request.POST.get("comment")
    mod_uid = request.POST.get("mod_uid")
    offtopic = request.POST.get("offtopic")
    is_valid, msg = False, "Invalid"

    if (action is None) and not (dupe or mod_uid):
        msg = "Select an action."
    if action == BUMP_POST and not post.is_toplevel:
        msg = "You can only perform this action to a top-level post"

    auth.moderate_post(post=post, request=request, action=action, comment=dupe_comment,
                       dupes=dupes, pid=mod_uid, offtopic=offtopic)

    parent = Post.objects.filter(uid=mod_uid).first()
    if not parent and mod_uid:
        msg = f"Parent id: {mod_uid} does not exist."

    if parent.root != post.root:
        msg = f"Parent does not share the same root."

    return is_valid, msg


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def ajax_moderate(request):
    user = request.user
    is_moderator = user.is_authenticated and user.profile.is_moderator
    uid = request.POST.get("uid")
    post = Post.objects.filter(uid=uid).first()

    # Bail out on errors
    if not is_moderator:
        return ajax_error(msg="You need to be a moderator to preform this action.")
    if not post:
        return ajax_error(msg=f"Post id: {uid} does not exist.")

    is_valid, msg = validate_moderation(request=request, post=post)

    return validate_moderation(request=request, post=post)


