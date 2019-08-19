from functools import wraps, partial
import logging
from ratelimit.decorators import ratelimit

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
def ajax_digest(request, uid):
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


@ajax_error_wrapper(method="GET")
def ajax_inplace(request):

    uid = request.GET.get("uid")
    post = Post.objects.filter(uid=uid).first()

    if not post:
        return ajax_error(msg="Post does not exist")
    rows = len(post.content.split("\n")) + 2
    tmpl = loader.get_template("widgets/inplace_form.html")
    # tmpl = loader.get_template("widgets/test_search_results.html")
    context = dict(post=post, rows=rows)

    inplace_form = tmpl.render(context)

    return ajax_success(msg="success", content=post.content, inplace_form=inplace_form)


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
