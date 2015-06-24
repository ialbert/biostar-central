from __future__ import absolute_import, division, print_function, unicode_literals

# Python modules.
import logging
from collections import OrderedDict
from django.db.models import Max, Count

# Django specific modules.
from django.shortcuts import render, redirect
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.db.models import Q, F
from django.contrib.auth.decorators import login_required
from haystack.query import SearchQuerySet
from django.core.cache import cache

from taggit.models import Tag

# Biostar specific local modules.
from . import query, search, auth
from .models import *
from biostar3.context import SESSION_COUNT_KEY

from biostar3.utils.compat import *

logger = logging.getLogger('biostar')

# Get custom user model.
User = get_user_model()


def tag_list(request):
    template_name = "tag_list.html"

    # Handle tag search parameter.
    q = request.GET.get('q', '')

    tags = Tag.objects
    tags = tags.filter(name__icontains=q) if q else tags.all()
    tags = tags.order_by("name")

    paginator = query.ExtendedPaginator(request, object_list=tags, time_class=query.DropDown,
                                        sort_class=query.TagSortValidator, per_page=100)
    page = paginator.curr_page()

    context = dict(page=page, tags=page.object_list, q=q)
    return render(request, template_name, context)


def tag_filter(request, name):
    """
    Returns a list of posts filtered by a tag name.
    """
    posts = query.get_toplevel_posts(user=request.user)
    names = name.split("+")
    posts = posts.filter(tags__name__in=names)
    messages.info(request, 'Filtering for tags: %s' % name)
    return post_list(request, posts=posts)


def unanswered(request):
    """
    Returns a list of posts filtered by a tag name.
    """
    posts = query.get_toplevel_posts(user=request.user).filter(type=Post.QUESTION, reply_count=0)
    messages.info(request, 'Unanswered questions')
    return post_list(request, posts=posts)


@auth.valid_user
def posts_by_user(request, pk, target=None):
    """
    Returns the posts created by a user.
    """
    posts = query.get_all_posts(request, target=target)
    messages.info(request, 'Posts by: %s' % target.name)
    return post_list(request, posts=posts)


@auth.valid_user
def upvoted_posts(request, pk, target=None):
    """
    Returns the upvoted posts created by a user.
    """
    posts = query.get_posts_by_vote(user=target, vote_types=[Vote.BOOKMARK, Vote.UP])
    messages.info(request, 'Upvoted posts created by %s' % target.name)
    return post_list(request, posts=posts)


@login_required
def my_bookmarks(request):
    """
    Returns the bookmarks by a user.
    """
    user = request.user
    posts = query.get_my_bookmarks(user=user)
    messages.info(request, 'Bookmarks for: %s' % user.name)
    return post_list(request, posts=posts)


@login_required
def my_messages(request):
    """
    Returns the messages for a user.
    """
    template_name = "message_list.html"
    user = request.user

    posts = Message.objects.filter(user=user).select_related("body").order_by('-date')

    paginator = query.ExtendedPaginator(request,
                                        object_list=posts, per_page=25)
    page = paginator.curr_page()

    # Force into list before the update takes place.
    context = dict(page=page, notes=list(page.object_list))

    # Set all messages to read.
    posts.update(unread=False)

    # Reset the value for for messages.
    counts = request.session.get(SESSION_COUNT_KEY, {})
    counts['mesg_count'] = 0  # Hardcoded key must match that set in contex.py
    request.session[SESSION_COUNT_KEY] = counts

    return render(request, template_name, context)


@login_required
@auth.valid_user
def vote_list(request, pk, target=None):
    template_name = "vote_list.html"

    votes = Vote.objects.filter(post__author=target).select_related("post", "author").order_by("-date")

    # Set all votes to seen.
    votes.update(unread=False)

    paginator = query.ExtendedPaginator(request,
                                        object_list=votes, per_page=50)
    page = paginator.curr_page()

    context = dict(page=page, votes=page.object_list, target=target)

    # Reset the value for for messages.
    counts = request.session.get(SESSION_COUNT_KEY, {})
    counts['new_vote_count'] = 0  # Hardcoded key must match that set in contex.py
    request.session[SESSION_COUNT_KEY] = counts

    return render(request, template_name, context)


def site_filter(request, posts):
    user = request.user

    # Authenticated users on the main site filter by subscription.
    if request.site.id == settings.SITE_ID and request.subs:
        return posts.filter(site_id__in=request.subs)

    # Default is to filter by site.
    return posts.filter(site=request.site)


def post_list(request, posts=None):
    template_name = "post_list.html"

    if posts is None:
        # The view is generic and could be called prefilled with posts.
        posts = query.get_toplevel_posts(user=request.user)

    # Filter posts by the site the user is accessing.
    posts = site_filter(request, posts=posts)

    # Add the paginator.
    paginator = query.ExtendedPaginator(request,
                                        sort_class=query.PostSortValidator,
                                        time_class=query.TimeLimitValidator,
                                        object_list=posts, per_page=settings.POSTS_PER_PAGE)

    page = paginator.curr_page()

    html_title = "Post List"
    context = dict(page=page, posts=page.object_list, html_title=html_title)

    return render(request, template_name, context)


def site_list(request):
    template_name = "site_list.html"

    user = request.user
    # Get the subscriptions for the user.
    if user.is_anonymous():
        subs = []
    else:
        subs = set([sub.site.id for sub in SiteSub.objects.filter(user=user)])

    if request.method == "POST":
        if user.is_anonymous():
            messages.error(request, "Please log in to use this feature!")
        else:
            # Handles site subscription
            site_ids = request.POST.getlist('site_id')
            # A simple sanity check
            site_ids = site_ids[:50]
            sites = Site.objects.filter(id__in=site_ids)
            SiteSub.objects.filter(user=user).delete()
            for site in sites:
                SiteSub.objects.create(user=user, site=site)

            subs = set([sub.site.id for sub in SiteSub.objects.filter(user=user)])
            # Delete the session if exists.
            del request.session[settings.SUBSCRIPTION_CACHE_NAME]

    sites = Site.objects.all().order_by("id")
    for site in sites:

        site.checked = site.id in subs

    html_title = "Site List"
    context = dict(sites=sites, html_title=html_title, site=request.site,
                   site_scheme=settings.SITE_SCHEME)

    return render(request, template_name, context)


def search_results(request):
    """
    Produces the search results
    """
    template_name = "post_search_results.html"
    q = request.GET.get('q', '')

    if not q:
        return redirect(reverse("home"))

    posts = search.plain(q)

    paginator = query.ExtendedPaginator(request, object_list=posts,
                                        per_page=settings.POSTS_PER_PAGE, orphans=False)

    page = paginator.curr_page()

    context = dict(page=page, posts=page.object_list, q=q)

    return render(request, template_name, context)


def update_post_views(request, post, minutes=settings.POST_VIEW_INTERVAL):
    """
    Views are updated per user session"
    """
    ip = auth.remote_ip(request)
    since = auth.ago(minutes=minutes)
    try:
        # One view per time interval from each IP address.
        if not PostView.objects.filter(ip=ip, post=post, date__gt=since):
            PostView.objects.create(ip=ip, post=post, date=right_now())
            Post.objects.filter(id=post.id).update(view_count=F('view_count') + 1)
    except Exception as exc:
        # Triggers if the IP address is spoofed and/or malformed.
        logger.error(exc)


@auth.post_view
def post_view(request, pk, post=None, user=None):
    """
    Generates the page that contains a full thread.
    """

    template_name = "post_detail.html"

    if not post.is_toplevel:
        # Post is not at top level. Redirect and and scroll the page to the right anchor.
        return redirect(post.get_absolute_url())

    # Update post views
    update_post_views(request=request, post=post)

    # Gets all objects in a thread. Moderators get deleted objects as well.
    thread = [px for px in query.get_thread(post, user)]

    # Collect votes for authenticated users
    store = {Vote.UP: set(), Vote.BOOKMARK: set()}

    if user.is_authenticated():
        # Authenticated users have votes only.
        pids = [p.id for p in thread]
        votes = Vote.objects.filter(post_id__in=pids, author=user).values_list("post_id", "type")
        for post_id, vote_type in votes:
            store.setdefault(vote_type, set()).add(post_id)

    # Extra attributes carry context into templates.
    post.upvotes = store[Vote.UP]
    post.bookmarks = store[Vote.BOOKMARK]

    # Set up additional attributes on each post
    write_access_check = auth.write_access_func(user=user)

    def decorator(p):
        # Each post needs to carry information on its status relative to the user.
        p.editable = write_access_check(post=p)
        p.has_vote = p.id in post.upvotes
        p.has_bookmark = p.id in post.bookmarks
        return p

    # Decorate all posts in the thread.
    thread = list(map(decorator, thread))

    # Decorate the main post as well.
    post = decorator(post)

    # Store answers in a separate list for simpler access.
    post.answers = filter(lambda p: p.type == Post.ANSWER, thread)
    post.answers = list(post.answers)

    # Comments will be stored in a dictionary for fast access.
    comment_list = filter(lambda pc: pc.type == Post.COMMENT, thread)

    # Collect comments into a dictionary keyed by the parent id with
    # comments as a value list
    post.comments = OrderedDict()
    for comment in comment_list:
        post.comments.setdefault(comment.parent.id, []).append(comment)

    # Get related objects
    related = SearchQuerySet().more_like_this(post)[:25]
    related = filter(lambda x: x.object and x.object.is_toplevel, related)

    # Add object to the context.
    context = dict(post=post, related=related)

    return render(request, template_name, context)


def planet_list(request):
    template_name = "planet_list.html"

    q = request.GET.get('q', '')

    if q:
        posts = search.plain(q, "forum.blogpost")
    else:
        posts = BlogPost.objects.order_by("creation_date")

    # Get the blogs in updated order.
    blogs = Blog.objects.all()
    blogs = blogs.annotate(updated_date=Max("blogpost__creation_date"),
                           count=Count("blogpost__id")).order_by("-updated_date", "-list_order")

    paginator = query.ExtendedPaginator(request,
                                        object_list=posts, per_page=settings.POSTS_PER_PAGE)
    page = paginator.curr_page()

    context = dict(posts=posts, page=page, blogs=blogs, q=q)

    return render(request, template_name, context)


from django.http import Http404


def flatpage_view(request, slug, domain=None, flatpage=None, user=None):
    domain = domain or settings.DEFAULT_GROUP_DOMAIN

    template_name = "page_view.html"

    flatpage = FlatPage.objects.filter(slug=slug).select_related("post", "author").first()

    if not flatpage:
        msg = "The page {}/{} does not seem to exists".format(domain, slug)
        raise Http404(msg)

    context = dict(flatpage=flatpage)
    return render(request, template_name, context)
