from functools import wraps, partial
import logging
import json
from ratelimit.decorators import ratelimit
from urllib import request as builtin_request
# import requests
from urllib.parse import urlencode
from datetime import datetime, timedelta
from urllib.parse import quote

from django.conf import settings
from django.db.models import Q, Count
from django.shortcuts import reverse, redirect
from django.template import loader
from django.http import JsonResponse
from django.utils.decorators import available_attrs
from django.db.models import F
from whoosh.searching import Results
from whoosh.sorting import FieldFacet, ScoreFacet
from .const import *
from taggit.models import Tag
from biostar.accounts.models import Profile, User
from . import auth, util, forms, tasks, search, views
from .models import Post, Vote, Subscription


def ajax_msg(msg, status, **kwargs):
    payload = dict(status=status, msg=msg)
    payload.update(kwargs)
    return JsonResponse(payload)


logger = logging.getLogger("biostar")
ajax_success = partial(ajax_msg, status='success')
ajax_error = partial(ajax_msg, status='error')

MIN_TITLE_CHARS = 10
MAX_TITLE_CHARS = 180

MAX_TAGS = 5


class ajax_error_wrapper:
    """
    Used as decorator to trap/display  errors in the ajax calls
    """

    def __init__(self, method, login_required=True):
        self.method = method
        self.login_required = login_required

    def __call__(self, func, *args, **kwargs):

        @wraps(func, assigned=available_attrs(func))
        def _ajax_view(request, *args, **kwargs):

            if request.method != self.method:
                return ajax_error(f'{self.method} method must be used.')

            if not request.user.is_authenticated and self.login_required:
                return ajax_error('You must be logged in.')

            return func(request, *args, **kwargs)

        return _ajax_view


def ajax_test(request):
    """
    Creates a commment on a top level post.
    """
    msg = "OK"
    print(f"HeRe= {request.POST} ")
    return ajax_error(msg=msg)


@ajax_error_wrapper(method="GET")
def user_image(request, username):
    user = User.objects.filter(username=username).first()

    gravatar_url = auth.gravatar(user=user)
    return redirect(gravatar_url)


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

    not_moderator = user.is_authenticated and not user.profile.is_moderator
    if post.root.author != user and not_moderator and vote_type == Vote.ACCEPT:
        return ajax_error("Only moderators or the person asking the question may accept answers.")

    msg, vote, change = auth.apply_vote(post=post, user=user, vote_type=vote_type)

    return ajax_success(msg=msg, change=change)


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST", login_required=True)
def drag_and_drop(request):
    was_limited = getattr(request, 'limited', False)
    if was_limited:
        return ajax_error(msg="Too many request from same IP address. Temporary ban.")

    parent_uid = request.POST.get("parent", '')
    uid = request.POST.get("uid", '')

    parent = Post.objects.filter(uid=parent_uid).first()
    post = Post.objects.filter(uid=uid).first()
    post_type = Post.COMMENT
    if not post:
        return ajax_error(msg="Post does not exist.")

    if parent_uid == "NEW":
        parent = post.root
        post_type = Post.ANSWER

    if not uid or not parent:
        return ajax_error(msg="Parent and Uid need to be provided. ")

    if not request.user.profile.is_moderator:
        return ajax_error(msg="Only moderators can move comments.")

    children = auth.walk_down_thread(parent=post)

    if parent == post or (parent in children):
        return ajax_error(msg="Can not move post under parent.")

    if post.is_toplevel:
        return ajax_error(msg="Top level posts can not be moved.")

    Post.objects.filter(uid=post.uid).update(type=post_type, parent=parent)

    post.update_parent_counts()
    redir = post.get_absolute_url()

    return ajax_success(msg="success", redir=redir)


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def ajax_subs(request):
    was_limited = getattr(request, 'limited', False)

    if was_limited:
        return ajax_error(msg="Too many request from same IP address. Temporary ban.")

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
def ajax_digest(request):
    was_limited = getattr(request, 'limited', False)
    user = request.user
    if was_limited:
        return ajax_error(msg="Too many request from same IP address. Temporary ban.")
    type_map = dict(daily=Profile.DAILY_DIGEST, weekly=Profile.WEEKLY_DIGEST,
                    monthly=Profile.MONTHLY_DIGEST)
    if user.is_anonymous:
        return ajax_error(msg="You need to be logged in to edit your profile.")

    # Get the post and digest preference.
    pref = request.POST.get('pref')
    pref = type_map.get(pref, Profile.NO_DIGEST)
    Profile.objects.filter(user=user).update(digest_prefs=pref)

    return ajax_success(msg="Changed digest options.")


def validate_recaptcha(token):
    """
    Send recaptcha token to API to check if user response is valid
    """
    url = 'https://www.google.com/recaptcha/api/siteverify'
    values = {
        'secret': settings.RECAPTCHA_PRIVATE_KEY,
        'response': token
    }
    data = urlencode(values).encode("utf-8")
    response = builtin_request.urlopen(url, data)
    result = json.load(response)

    if result['success']:
        return True, ""

    return False, "Invalid reCAPTCHA. Please try again."


def validate_toplevel(title, post_type, tags_list, parent=None, user=None):
    """Validate fields found in top level posts"""

    title_length = len(title.replace(' ', ''))
    allowed_types = [opt[0] for opt in Post.TYPE_CHOICES]
    tag_length = len(tags_list)

    if title_length <= MIN_TITLE_CHARS:
        msg = f"Title too short, please add more than {MIN_TITLE_CHARS} characters."
        return False, msg
    if title_length > MAX_TITLE_CHARS:
        msg = f"Title too long, please add less than {MAX_TITLE_CHARS} characters."
        return False, msg

    if post_type not in allowed_types:
        msg = "Not a valid post type."
        return False, msg

    if tag_length > MAX_TAGS:
        msg = f"Too many tags, maximum of {MAX_TAGS} tags allowed."
        return False, msg

    if parent and user:
        if parent.root.is_locked and not user.is_superuser:
            msg = "This post is locked. Only admins can contribute to it."
            return False, msg

        if parent.root.closed and (user.is_anonymous or not user.profile.is_moderator or parent.root.owner != user):
            msg = "This post is closed. Only moderators and the initial author can contribute to it."
            return False, msg

    return True, ""


def validate_post(content, title, tags_list, post_type, is_toplevel=False, recaptcha_token='',
                  check_captcha=False,  parent=None, user=None):
    content_length = len(content.replace(' ', ''))

    if check_captcha:
        valid_captcha, msg = validate_recaptcha(recaptcha_token)
        if not valid_captcha:
            return False, msg

    # Validate fields common to all posts.
    if content_length <= forms.MIN_CONTENT:
        msg = f"Content too short, please add more than add more {forms.MIN_CONTENT} characters."
        return False, msg
    if content_length > forms.MAX_CONTENT:
        msg = f"Content too long, please add less than {forms.MAX_CONTENT} characters."
        return False, msg

    # Validate fields found in top level posts
    if is_toplevel:
        return validate_toplevel(title=title, post_type=post_type, tags_list=tags_list, parent=parent,
                                 user=user)

    return True, ""


def is_trusted(user):
    # Moderators and users with scores above threshold are trusted.
    trusted = user.is_authenticated and (user.profile.trusted or user.profile.score > 15)
    return trusted


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def ajax_edit(request, uid):
    """
    Edit post content using ajax.
    """
    was_limited = getattr(request, 'limited', False)
    if was_limited:
        return ajax_error(msg="Too many request from same IP address. Temporary ban.")

    post = Post.objects.filter(uid=uid).first()
    if not post:
        return ajax_error(msg="Post does not exist")

    content = request.POST.get("content", post.content)
    title = request.POST.get("title", post.title)
    post_type = int(request.POST.get("type", post.type))
    tag_list = set(request.POST.getlist("tag_val", []))
    tag_str = ','.join(tag_list)

    # Validate fields in request.POST
    valid, msg = validate_post(content=content, title=title, tags_list=tag_list,
                               post_type=post_type, is_toplevel=post.is_toplevel)
    if not valid:
        return ajax_error(msg=msg)

    if post.is_toplevel:
        post.title = title
        post.type = post_type
        post.tag_val = tag_str

    post.lastedit_user = request.user
    post.content = content
    post.save()

    tags = post.tag_val.split(",")
    context = dict(post=post, tags=tags, show_views=True)

    tmpl = loader.get_template('widgets/post_tags.html')
    tag_html = tmpl.render(context)
    new_title = f'{post.get_type_display()}: {post.title}'

    return ajax_success(msg='success', html=post.html, title=new_title, tag_html=tag_html)


# @ajax_error_wrapper(method="GET")
# def ajax_recent(request):
#
#     uid = request.GET.get('uid', '')
#     # Return recent posts made in past 5000 seconds.
#     seconds = 500000
#     delta = util.now() - timedelta(seconds=seconds)
#     if uid:
#         root = Post.objects.filter(uid=uid).first().root
#         new_posts = Post.objects.filter(root=root, lastedit_date__gt=delta).exclude(uid=uid)
#     else:
#         new_posts = Post.objects.filter(type__in=Post.TOP_LEVEL, lastedit_date__gt=delta)
#
#     if not new_posts:
#         return ajax_success(msg='No new posts in the past 5000 seconds.', data='')
#
#     nposts = len(new_posts)
#     context = dict(nposts=nposts)
#     tmpl = loader.get_template('widgets/ajax_recent.html')
#     template = tmpl.render(context=context)
#     most_recent = new_posts.order_by('-pk').first()
#     most_recent_url = most_recent.get_absolute_url() if most_recent is not None else ''
#
#     return ajax_success(msg="New posts", template=template)


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def ajax_create(request):
    was_limited = getattr(request, 'limited', False)
    if was_limited:
        return ajax_error(msg="Too many request from same IP address. Temporary ban.")

    # Get form fields from POST request
    user = request.user
    content = request.POST.get("content", '')
    title = request.POST.get("title", '')
    tag_list = {x.strip() for x in request.POST.getlist("tag_val", [])}

    tag_str = ','.join(tag_list)
    recaptcha_token = request.POST.get("recaptcha_response")
    is_toplevel = bool(int(request.POST.get('top', 0)))
    # Get the post type
    post_type = request.POST.get('type', '0')
    post_type = int(post_type) if post_type.isdigit() else 0

    # Find out if we are currently creating a comment
    is_comment = request.POST.get('comment', '0')
    is_comment = int(is_comment) if is_comment.isdigit() else 0

    # Resolve the parent post
    parent_uid = request.POST.get("parent", '')
    parent = Post.objects.filter(uid=parent_uid).first()

    # reCAPTCHA field required when an untrusted user creates any post.
    check_captcha = not is_trusted(user=user) and settings.RECAPTCHA_PRIVATE_KEY

    # Validate fields in request.POST
    valid, msg = validate_post(content=content, title=title, recaptcha_token=recaptcha_token, tags_list=tag_list,
                               post_type=post_type, check_captcha=check_captcha, is_toplevel=is_toplevel)

    if not valid:
        return ajax_error(msg=msg)

    if parent_uid and not parent:
        return ajax_error(msg='Parent post does not exist.')

    # We are creating an answer.
    if parent and parent.is_toplevel and not is_comment:
        post_type = Post.ANSWER
    # Creating a comment
    elif is_comment:
        post_type = Post.COMMENT

    # Create the post.
    post = Post.objects.create(title=title, tag_val=tag_str, type=post_type, content=content,
                               author=user, parent=parent)
    return ajax_success(msg='Created post', redirect=post.get_absolute_url())


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="GET")
def inplace_form(request):
    user = request.user
    if user.is_anonymous:
        return ajax_error(msg="You need to be logged in to edit or create posts.")

    # Get the current post uid we are editing
    uid = request.GET.get('uid', '')
    # Parent uid to add an answer/comment to.
    parent_uid = request.GET.get('parent', '')
    post = Post.objects.filter(uid=uid).first()
    if uid and not post:
        return ajax_error(msg="Post does not exist.")

    is_toplevel = post.is_toplevel if post else int(request.GET.get('top', 1))
    is_comment = request.GET.get('comment', 0) if not post else post.is_comment

    title = post.title if post else ''
    content = post.content if post else ''
    rows = len(post.content.split("\n")) if post else request.GET.get('rows', 10)
    rows = rows if 25 > int(rows) >= 2 else 4
    html = post.html if post else ''

    template = "widgets/inplace_form.html"

    # Untrusted users get a reCAPTCHA field added when creating posts.
    creating = post is None
    not_trusted = not is_trusted(user=user) if creating else False

    # Load the content and form template
    tmpl = loader.get_template(template_name=template)
    context = dict(user=user, post=post, content=content, rows=rows, is_comment=is_comment, is_toplevel=is_toplevel,
                   parent_uid=parent_uid, captcha_key=settings.RECAPTCHA_PUBLIC_KEY, html=html,
                   not_trusted=not_trusted, title=title)

    form = tmpl.render(context)
    return ajax_success(msg="success", inplace_form=form)


def close(r):
    # Ensure the searcher object gets closed.
    r.searcher.close() if isinstance(r, Results) else None
    return


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
def ajax_search(request):
    query = request.GET.get('query', '')
    page = request.GET.get('page', 1)

    try:
        redir = bool(int(request.GET.get('redir', 0)))
    except Exception as exc:
        redir = False

    if redir:
        # Redirect search results to a separate page.
        redir_url = reverse('post_search') + '?query=' + query
        return ajax_success(redir=redir_url, msg="success")

    if not query:
        return ajax_success(msg="Empty query", status="error")

    # Preform search on indexed posts.
    results = search.preform_whoosh_search(query=query, page=page, per_page=settings.SEARCH_RESULTS_PER_PAGE)

    tmpl = loader.get_template("widgets/search_results.html")
    context = dict(results=results, query=query)

    results_html = tmpl.render(context)

    return ajax_success(html=results_html, msg="success")


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
def ajax_tags_search(request):
    was_limited = getattr(request, 'limited', False)
    if was_limited:
        return ajax_error(msg="Too many request from same IP address. Temporary ban.")

    query = request.GET.get('query', '')
    count = Count('post', filter=Q(post__type__in=Post.TOP_LEVEL))

    if len(query) < settings.SEARCH_CHAR_MIN:
        return ajax_error(msg=f"Enter more than {settings.SEARCH_CHAR_MIN} characters")

    if query:
        db_query = Q(name__in=query) | Q(name__contains=query)

        results = Tag.objects.annotate(tagged=count).order_by('-tagged').filter(db_query)

        tmpl = loader.get_template("widgets/search_results.html")
        context = dict(results=results, query=query, tags=True)
        results_html = tmpl.render(context)
        return ajax_success(html=results_html, msg="success")

    return ajax_success(msg="success")


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
def ajax_users_search(request):
    was_limited = getattr(request, 'limited', False)
    if was_limited:
        return ajax_error(msg="Too many request from same IP address. Temporary ban.")

    query = request.GET.get('query', '')
    if len(query) < settings.SEARCH_CHAR_MIN:
        return ajax_error(msg=f"Enter more than {settings.SEARCH_CHAR_MIN} characters")
    redir = reverse('community_list')
    redir = redir + "?query=" + quote(query)

    return ajax_success(redir=redir, msg="success")


def similar_posts(request, uid):
    """
    Return a feed populated with posts similar to the one in the request.
    """

    post = Post.objects.filter(uid=uid).first()
    if not post:
        ajax_error(msg='Post does not exist.')

    results = search.more_like_this(uid=post.uid)

    tmpl = loader.get_template('widgets/similar_posts.html')
    context = dict(results=results)
    results_html = tmpl.render(context)

    return ajax_success(html=results_html, msg="success")
