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


def add_to_closed(post, user):

    if user.is_moderator or post.author == user:
        return True


@ajax_error_wrapper(method="GET", login_required=True)
def report_spammer(request, post_uid):
    """
    Report this user as a spammer.
    """

    post = Post.objects.filter(uid=post_uid).first()

    if not post:
        return ajax_error(msg='Post does not exist.')

    if request.user == post.author or post.author.profile.is_moderator:
        return ajax_error(msg='Invalid action.')

    auth.handle_spam_post(post=post)

    return ajax_success(msg="Reported user as a spammer.")


def validate_root(post, user):

    if not post:
        return False, "No post provided."

    # Anonymous users and regular users who are not the author
    # are not allowed to add to closed posts.
    allowed = user.is_authenticated and (user.profile.is_moderator or post.root.author == user)
    if not post.root.is_open and not allowed:
        msg = "This post is not open. Only moderators and the initial author can contribute to it."
        return False, msg

    return True, ""


def validate_toplevel(fields={}):
    """Validate fields found in top level posts"""

    title = fields.get('title', '')
    tags_list = fields.get('tags_list', '')
    post_type = fields.get('post_type', '')
    user = fields.get('user', '')
    root = fields.get('parent', '')

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

    if root and user:
        return validate_root(post=root, user=user)

    return True, ""


def validate_post(fields={}):
    """
    Validate fields found in dictionary.
    """
    content = fields.get('content', '')
    content_length = len(content.replace(' ', ''))
    recaptcha_token = fields.get('recaptcha_token')

    if fields.get('check_captcha'):
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
    if fields.get('is_toplevel'):
        return validate_toplevel(fields=fields)

    return True, ""


def moderate_user(request, target_id):

    user = request.user

    target_uid = User.objects.filter(id=target_id).first()

    return


def is_trusted(user):
    # Moderators and users with scores above threshold are trusted.
    trusted = user.is_authenticated and (user.profile.trusted or
                                         user.profile.score > settings.RECAPTCHA_TRUSTED_USER_SCORE)
    return trusted


def get_fields(request, post=None):
    """
    Used to retrieve all fields in a request used for editing and creating posts.
    """

    def request_data(data_key, placeholder=None):
        return request.POST.get(data_key, request.GET.get(data_key, placeholder))

    # Fields common to all posts
    user = request.user
    content = request_data("content", "")
    title = request_data("title", "")
    post_type = request_data("type", Post.QUESTION)
    post_type = int(post_type) if str(post_type).isdigit() else Post.QUESTION
    root = None

    if post:
        content = content or post.content
        title = title or post.title
        post_type = post_type if post_type is not None else Post.QUESTION
        root = post.root

    parent_uid = request_data('parent', '')
    parent = Post.objects.filter(uid=parent_uid).first()
    recaptcha_token = request_data("recaptcha_response")
    # reCAPTCHA field required when an untrusted user creates any post.
    check_captcha = not is_trusted(user=user) and settings.RECAPTCHA_PRIVATE_KEY != ''

    # Fields found in top level posts
    tag_list = {x.strip() for x in request.POST.getlist("tag_val", [])}
    tag_str = ','.join(tag_list)

    # Field flags. TODO: being refactored out.
    is_toplevel = bool(int(request_data('top', 0)))
    is_comment = request_data('comment', '0')
    is_comment = int(is_comment) if is_comment.isdigit() else 0

    fields = dict(content=content, title=title, post_type=post_type, tag_list=tag_list, tag_str=tag_str,
                  user=user, root=root, parent=parent, parent_uid=parent_uid, check_captcha=check_captcha,
                  recaptcha_token=recaptcha_token, is_comment=is_comment, is_toplevel=is_toplevel)

    return fields


def set_post(fields, post, save=True):

    if post.is_toplevel:
        post.title = fields.get('title', post.title)
        post.type = fields.get('post_type', post.type)
        post.tag_val = fields.get('tag_str', post.tag_val)

    post.lastedit_user = fields.get('user', post.lastedit_user)
    post.content = fields.get('content', post.content)
    if save:
        post.save()

    return post


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

    # Get the fields found in the request
    fields = get_fields(request=request, post=post)
    # Validate fields in request.POST
    valid, msg = validate_post(fields=fields)
    if not valid:
        return ajax_error(msg=msg)

    # Set the fields for this post.
    post = set_post(post=post, fields=fields)

    # Get the new tags and render them
    tags = post.tag_val.split(",")
    context = dict(post=post, tags=tags, show_views=True)
    tmpl = loader.get_template('widgets/post_tags.html')
    tag_html = tmpl.render(context)

    # Prepare the new title
    new_title = f'{post.get_type_display()}: {post.title}'

    return ajax_success(msg='success', html=post.html, title=new_title, tag_html=tag_html)


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def ajax_create(request):
    was_limited = getattr(request, 'limited', False)
    if was_limited:
        return ajax_error(msg="Too many request from same IP address. Temporary ban.")

    # Get form fields from POST request
    fields = get_fields(request=request, post=None)
    user = request.user

    # Validate fields in request.POST
    valid, msg = validate_post(fields=fields)

    if not valid:
        return ajax_error(msg=msg)
    parent = fields.get('parent')
    parent_uid = fields.get('parent_uid')
    is_comment = fields.get('is_comment')

    if parent_uid and not parent:
        return ajax_error(msg='Parent post does not exist.')

    # We are creating an answer.
    if parent and parent.is_toplevel and not is_comment:
        post_type = Post.ANSWER
    # Creating a comment
    elif is_comment:
        post_type = Post.COMMENT
    else:
        post_type = fields.get('post_type')

    # Create the post.
    post = Post.objects.create(title=fields.get('title'), tag_val=fields.get('tag_str'), type=post_type,
                               content=fields.get('content'), author=user, parent=fields.get('parent'))

    return ajax_success(msg='Created post', redirect=post.get_absolute_url())


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="GET")
def inplace_form(request):
    """
    Used to render inplace forms for creation and editting posts.
    """
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

    # Extract the
    title = post.title if post else ''
    content = post.content if post else ''
    rows = len(post.content.split("\n")) if post else request.GET.get('rows', 10)
    rows = rows if 25 > int(rows) >= 2 else 4
    html = post.html if post else ''

    template = "widgets/inplace_form.html"

    # Untrusted users get a reCAPTCHA field added when creating posts.
    creating = post is None
    parent = Post.objects.filter(uid=parent_uid).first()
    not_trusted = not is_trusted(user=user) if creating else False
    if parent or post:
        # Validate that the root is not locked or closed.
        valid, msg = validate_root(post=parent or post, user=user)
        if not valid:
            return ajax_error(msg=msg)

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

    # if settings.USE_ELASTIC_SEARCH:
    #     results = search.preform_elastic_search(post=post)
    #     template_name = "widgets/similar_posts_elastic.html"
    #     print(results)
    #
    # else:
    results = search.whoosh_more_like_this(uid=post.uid)
    template_name = 'widgets/similar_posts.html'

    tmpl = loader.get_template(template_name)
    context = dict(results=results)
    results_html = tmpl.render(context)

    return ajax_success(html=results_html, msg="success")
