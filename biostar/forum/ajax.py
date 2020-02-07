from functools import wraps, partial
import logging
import json
from ratelimit.decorators import ratelimit
from urllib import request as builtin_request
# import requests
from urllib.parse import urlencode
from datetime import datetime, timedelta
from urllib.parse import quote

from django.core.cache import cache
from django.conf import settings
from django.db.models import Q, Count
from django.shortcuts import reverse, redirect
from django.template import loader
from django.http import JsonResponse

from whoosh.searching import Results

from biostar.accounts.models import Profile, User
from . import auth, util, forms, tasks, search, views, const
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


class ajax_error_wrapper:
    """
    Used as decorator to trap/display  errors in the ajax calls
    """

    def __init__(self, method, login_required=True):
        self.method = method
        self.login_required = login_required

    def __call__(self, func, *args, **kwargs):

        @wraps(func)
        def _ajax_view(request, *args, **kwargs):

            if request.method != self.method:
                return ajax_error(f'{self.method} method must be used.')

            if not request.user.is_authenticated and self.login_required:
                return ajax_error('You must be logged in.')

            if request.user.is_authenticated and request.user.profile.is_spammer:
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
    user = request.user
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

    if not (user.profile.is_moderator or post.author == user):
        return ajax_error(msg="Only moderators or the author can move posts.")

    collect = set()
    auth.walk_down_thread(parent=post, collect=collect)

    if parent == post or (parent in collect) or parent.root != post.root:
        return ajax_error(msg="Can not move post here.")

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


@ajax_error_wrapper(method="GET", login_required=True)
def most_recent_users(request):

    recent_locations = Profile.objects.exclude(Q(location="") | Q(state__in=[Profile.BANNED, Profile.SUSPENDED])).prefetch_related("user")
    recent_locations = recent_locations.order_by('last_login')
    recent_locations = recent_locations[:settings.LOCATION_FEED_COUNT]

    users = [dict(username=u.user.username, email=u.user.email, uid=u.uid, name=u.name,
                  url=u.get_absolute_url(), score=u.score,
                  gravatar=auth.gravatar(user=u.user, size=30))
             for u in recent_locations]

    return ajax_success(msg="recent users", users=users)


@ajax_error_wrapper(method="GET", login_required=True)
def report_spam(request, post_uid):
    """
    Report this user as a spammer.
    """

    post = Post.objects.filter(uid=post_uid).first()

    if not post:
        return ajax_error(msg='Post does not exist.')

    if not request.user.profile.is_moderator:
        return ajax_error(msg="You need to be a moderator to preform that action.")

    if request.user == post.author or post.author.profile.is_moderator:
        return ajax_error(msg='Invalid action.')

    auth.handle_spam_post(post=post, user=request.user)

    return ajax_success(msg="Reported user as a spammer.")


def validate_toplevel_fields(fields={}):
    """Validate fields found in top level posts"""

    title = fields.get('title', '')
    tag_list = fields.get('tag_list', [])
    tag_val = fields.get('tag_val', '')
    post_type = fields.get('post_type', '')
    title_length = len(title.replace(' ', ''))
    allowed_types = [opt[0] for opt in Post.TYPE_CHOICES]
    tag_length = len(tag_list)

    if title_length <= forms.MIN_CHARS:
        msg = f"Title too short, please add more than {forms.MIN_CHARS} characters."
        return False, msg

    if title_length > forms.MAX_TITLE:
        msg = f"Title too long, please add less than {forms.MAX_TITLE} characters."
        return False, msg

    if post_type not in allowed_types:
        msg = "Not a valid post type."
        return False, msg
    if tag_length > forms.MAX_TAGS:
        msg = f"Too many tags, maximum of {forms.MAX_TAGS} tags allowed."
        return False, msg

    if len(tag_val) > 100:
        msg = f"Tags have too many characters, maximum of 100 characters total allowed in tags."
        return False, msg

    return True, ""


def validate_post_fields(fields={}, is_toplevel=False):
    """
    Validate fields found in dictionary.
    """
    content = fields.get('content', '')
    content_length = len(content.replace(' ', ''))
    user = fields.get("user")
    recaptcha_token = fields.get('recaptcha_token')
    # reCAPTCHA field required when an untrusted user creates any post.
    check_captcha = user.profile.require_recaptcha() and settings.RECAPTCHA_PRIVATE_KEY != ''

    if check_captcha:
        valid_captcha, msg = validate_recaptcha(recaptcha_token)
        if not valid_captcha:
            return False, msg

    # Validate fields common to all posts.
    if content_length <= forms.MIN_CHARS:
        msg = f"Content too short, please add more than {forms.MIN_CHARS} characters."
        return False, msg
    if content_length > forms.MAX_CONTENT:
        msg = f"Content too long, please add less than {forms.MAX_CONTENT} characters."
        return False, msg

    # Validate fields found in top level posts
    if is_toplevel:
        return validate_toplevel_fields(fields=fields)

    return True, ""


def get_fields(request, post=None):
    """
    Used to retrieve all fields in a request used for editing and creating posts.
    """
    user = request.user
    content = request.POST.get("content", "") or post.content
    title = request.POST.get("title", "") or post.title
    post_type = request.POST.get("type", Post.QUESTION)

    post_type = int(post_type) if str(post_type).isdigit() else Post.QUESTION
    recaptcha_token = request.POST.get("recaptcha_response")

    # Fields found in top level posts
    tag_list = {x.strip() for x in request.POST.getlist("tag_val", [])}
    tag_val = ','.join(tag_list) or post.tag_val

    fields = dict(content=content, title=title, post_type=post_type, tag_list=tag_list, tag_val=tag_val,
                  user=user, recaptcha_token=recaptcha_token,)

    return fields


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST", login_required=True)
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

    # Get the fields found in the request
    fields = get_fields(request=request, post=post)

    if not (request.user.profile.is_moderator or request.user == post.author):
        return ajax_error(msg="Only moderators or the author can edit posts.")

    # Validate fields in request.POST
    valid, msg = validate_post_fields(fields=fields, is_toplevel=post.is_toplevel)
    if not valid:
        return ajax_error(msg=msg)

    # Set the fields for this post.
    if post.is_toplevel:
        post.title = fields.get('title', post.title)
        post.type = fields.get('post_type', post.type)
        post.tag_val = fields.get('tag_val', post.tag_val)
    post.lastedit_user = request.user
    post.lastedit_date = util.now()
    post.content = fields.get('content', post.content)
    post.save()

    # Get the newly set tags to render
    tags = post.tag_val.split(",")
    context = dict(post=post, tags=tags, show_views=True)
    tmpl = loader.get_template('widgets/post_tags.html')
    tag_html = tmpl.render(context)

    # Prepare the new title to render
    new_title = f'{post.get_type_display()}: {post.title}'

    return ajax_success(msg='success', html=post.html, title=new_title, tag_html=tag_html)


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def ajax_comment_create(request):
    was_limited = getattr(request, 'limited', False)
    if was_limited:
        return ajax_error(msg="Too many request from same IP address. Temporary ban.")

    # Fields common to all posts
    user = request.user
    content = request.POST.get("content", "")
    recaptcha_token = request.POST.get("recaptcha_response")

    parent_uid = request.POST.get('parent', '')
    parent = Post.objects.filter(uid=parent_uid).first()
    if not parent:
        return ajax_error(msg='Parent post does not exist.')

    fields = dict(content=content, user=user, parent=parent, recaptcha_token=recaptcha_token)
    # Validate the fields
    valid, msg = validate_post_fields(fields=fields, is_toplevel=False)
    if not valid:
        return ajax_error(msg=msg)

    # Create the comment.
    post = Post.objects.create(type=Post.COMMENT, content=content, author=user, parent=parent)

    return ajax_success(msg='Created post', redirect=post.get_absolute_url())


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="GET")
def inplace_form(request):
    """
    Used to render inplace forms for editing posts.
    """
    user = request.user
    if user.is_anonymous:
        return ajax_error(msg="You need to be logged in to edit or create posts.")

    # Get the current post uid we are editing
    uid = request.GET.get('uid', '')
    post = Post.objects.filter(uid=uid).first()
    if not post:
        return ajax_error(msg="Post does not exist.")

    add_comment = request.GET.get('add_comment', False)

    # Load the content and form template
    template = "forms/inplace_form.html"
    tmpl = loader.get_template(template_name=template)
    users_str = auth.get_users_str()
    context = dict(user=user, post=post, rows=25, add_comment=add_comment, users_str=users_str,
                   captcha_key=settings.RECAPTCHA_PUBLIC_KEY)
    form = tmpl.render(context)

    return ajax_success(msg="success", inplace_form=form)


def similar_posts(request, uid):
    """
    Return a feed populated with posts similar to the one in the request.
    """

    post = Post.objects.filter(uid=uid).first()
    if not post:
        ajax_error(msg='Post does not exist.')

    cache_key = f"{const.SIMILAR_CACHE_KEY}-{post.uid}"
    results = cache.get(cache_key)

    if results is None:
        logger.info("Setting similar posts cache.")

        results = search.whoosh_more_like_this(uid=post.uid)
        # Set the results cache for 1 hour
        cache.set(cache_key, results, 3600)

    template_name = 'widgets/similar_posts.html'

    tmpl = loader.get_template(template_name)
    context = dict(results=results)
    results_html = tmpl.render(context)

    return ajax_success(html=results_html, msg="success")
