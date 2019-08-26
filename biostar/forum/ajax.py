from functools import wraps, partial
import logging
import json
from ratelimit.decorators import ratelimit
from urllib import request as builtin_request
#import requests
from urllib.parse import urlencode

from django.conf import settings
from django.db.models import Q, Count
from django.template import loader
from django.http import JsonResponse
from django.utils.decorators import available_attrs
from whoosh.searching import Results
from whoosh.sorting import FieldFacet, ScoreFacet
from .const import *
from taggit.models import Tag
from biostar.accounts.models import Profile
from . import auth, util, forms, tasks, search
from .models import Post, Vote, Subscription


def ajax_msg(msg, status, **kwargs):
    payload = dict(status=status, msg=msg)
    payload.update(kwargs)
    return JsonResponse(payload)


logger = logging.getLogger("biostar")
ajax_success = partial(ajax_msg, status='success')
ajax_error = partial(ajax_msg, status='error')


MIN_TITLE_CHARS = 10
MAX_TITLE_CHARS = 5000

MAX_TAGS = 5

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

    not_moderator = user.is_authenticated and not user.profile.is_moderator
    if post.root.author != user and not_moderator and vote_type == Vote.ACCEPT:
        return ajax_error("Only moderators or the person asking the question may accept answers.")

    msg, vote, change = auth.apply_vote(post=post, user=user, vote_type=vote_type)

    return ajax_success(msg=msg, change=change)


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


def validate_post(content, title, tags_list, post_type, is_toplevel=False, recaptcha_token=None):
    content_length = len(content.replace(" ", ''))
    title_length = len(title.replace(' ', ''))
    allowed_types = [opt[0] for opt in Post.TYPE_CHOICES]
    tag_length = len(tags_list)

    if recaptcha_token:
        return validate_recaptcha(recaptcha_token)

    # Validate fields common to all posts.
    if content_length <= forms.MIN_CONTENT:
        msg = f"Content too short, please add more than add more {forms.MIN_CONTENT} characters."
        return False, msg
    if content_length > forms.MAX_CONTENT:
        msg = f"Content too long, please add less than {forms.MAX_CONTENT} characters."
        return False, msg

    # Validate fields found in top level posts
    if is_toplevel:
        if title_length <= MIN_TITLE_CHARS:
            msg = f"Title too short, please add more than add more {MIN_TITLE_CHARS} characters."
            return False, msg
        if title_length > MAX_TITLE_CHARS:
            msg = f"Title too long, please add more than add more {MAX_TITLE_CHARS} characters."
            return False, msg

        if post_type not in allowed_types:
            msg = "Not a valid post type."
            return False, msg

        if tag_length > MAX_TAGS:
            msg = f"Too many tags, maximum of {MAX_TAGS} tags allowed."
            return False, msg

    return True, ""


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
    tag_list = set(request.POST.getlist("tag_val", []))
    tag_str = ','.join(tag_list)
    recaptcha_token = request.POST.get("recaptcha_response")

    # Get the post type
    post_type = request.POST.get('type', '0')
    post_type = int(post_type) if post_type.isdigit() else 0

    # Find out if we are currently creating a comment
    is_comment = request.POST.get('comment', '0')
    is_comment = int(is_comment) if is_comment.isdigit() else 0

    # Resolve the parent post
    parent_uid = request.POST.get("parent", '')
    parent = Post.objects.filter(uid=parent_uid).first()

    # Validate fields in request.POST
    valid, msg = validate_post(content=content, title=title, recaptcha_token=recaptcha_token, tags_list=tag_list,
                               post_type=post_type)
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
    if is_toplevel:
        template = "widgets/edit_toplevel.html"
    else:
        template = "widgets/edit.html"

    title = post.title if post else ''
    content = post.content if post else ''
    rows = len(post.content.split("\n")) if post else request.GET.get('rows', 10)
    is_comment = request.GET.get('comment', 0)
    html = post.html if post else ''

    # Untrusted users get a reCAPTCHA field added to the form
    not_trusted = user.is_authenticated and (not user.profile.trusted)
    # Load the content and form template
    tmpl = loader.get_template(template_name=template)
    context = dict(user=user, post=post, content=content, rows=rows, is_comment=is_comment,
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
    fields = ['content', 'tag_val', 'title', 'author', 'author_uid', 'author_handle']

    if query:

        whoosh_results = search.search(query=query, fields=fields)
        results = sorted(whoosh_results, key=lambda x: x['lastedit_date'], reverse=True)

        tmpl = loader.get_template("widgets/search_results.html")
        #tmpl = loader.get_template("widgets/test_search_results.html")
        context = dict(results=results, query=query)

        results_html = tmpl.render(context)
        # Ensure the whoosh reader is closed
        close(whoosh_results)
        return ajax_success(html=results_html, msg="success")

    return ajax_success(msg="Empty query, Enter atleast", status="error")


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
def ajax_tags_search(request):

    query = request.GET.get('query', '')
    #fields = ['content', 'tag_val', 'title', 'author', 'author_uid', 'author_handle']
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

    return ajax_success(msg="")


def similar_posts(request, uid):
    """
    Return a feed populated with posts similar to the one in the request.
    """

    post = Post.objects.filter(uid=uid).first()
    if not post:
        ajax_error(msg='Post does not exist.')

    results = []
    # Retrieve this post from the search index.
    indexed_post = search.preform_search(query=post.uid, fields=['uid'])

    if isinstance(indexed_post, Results) and not indexed_post.is_empty():

        results = indexed_post[0].more_like_this("content", top=settings.SIMILAR_FEED_COUNT)
        # Filter results for toplevel posts.
        results = filter(lambda p: p['is_toplevel'] is True, results)
        # Sort by lastedit_date
        results = sorted(results, key=lambda x: x['lastedit_date'], reverse=True)

    tmpl = loader.get_template('widgets/similar_posts.html')
    context = dict(results=results)
    results_html = tmpl.render(context)
    # Ensure the searcher object gets closed.
    close(indexed_post)

    return ajax_success(html=results_html, msg="success")
