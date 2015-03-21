from __future__ import absolute_import, division, print_function, unicode_literals

# Python modules.
import logging
from collections import OrderedDict

# Django specific modules.
from django.shortcuts import render, redirect
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.db.models import Q, F
from django.contrib.auth.decorators import login_required
from haystack.query import SearchQuerySet

from taggit.models import Tag

# Biostar specific local modules.
from . import models, query, search, auth
from .models import Vote, Post, PostView, UserGroup, GroupSub, Message, GroupPerm
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
    posts = query.get_toplevel_posts(user=request.user, group=request.group)
    names = name.split("+")
    posts = posts.filter(tags__name__in=names)
    messages.info(request, 'Filtering for tags: %s' % name)
    return post_list(request, posts=posts)


def unanswered(request):
    """
    Returns a list of posts filtered by a tag name.
    """
    posts = query.get_toplevel_posts(user=request.user, group=request.group).filter(type=Post.QUESTION, reply_count=0)
    messages.info(request, 'Unanswered questions')
    return post_list(request, posts=posts)


@auth.valid_user
def posts_by_user(request, pk, target=None):
    """
    Returns the posts created by a user.
    """
    posts = query.get_all_posts(request, target=target, group=request.group)
    messages.info(request, 'Posts by: %s' % target.name)
    return post_list(request, posts=posts)


@auth.valid_user
def upvoted_posts(request, pk, target=None):
    """
    Returns the upvoted posts created by a user.
    """
    posts = query.get_posts_by_vote(user=target, group=request.group, vote_types=[Vote.BOOKMARK, Vote.UP])
    messages.info(request, 'Upvoted posts created by %s' % target.name)
    return post_list(request, posts=posts)


@login_required
def my_bookmarks(request):
    """
    Returns the bookmarks by a user.
    """
    user = request.user
    posts = query.get_my_bookmarks(user=user, group=request.group)
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


def post_list(request, posts=None):
    template_name = "post_list.html"

    if posts is None:
        # The view is generic and could be called prefilled with posts.
        posts = query.get_toplevel_posts(user=request.user, group=request.group)

    paginator = query.ExtendedPaginator(request,
                                        sort_class=query.PostSortValidator,
                                        time_class=query.TimeLimitValidator,
                                        object_list=posts, per_page=settings.POSTS_PER_PAGE)
    page = paginator.curr_page()

    # Add the recent votes
    recent_votes = query.recent_votes()
    html_title = "Post List"
    context = dict(page=page, posts=page.object_list, recent_votes=recent_votes, html_title=html_title)

    return render(request, template_name, context)


@auth.group_access
def group_login(request, pk, group=None, user=None):
    # Handler fired on signup redirect when logging in from subdomain.
    # It auto adds user to the group that initiated the login process.
    return group_redirect_handler(request=request, group=group, user=user, autoadd=True)


@auth.group_access
def group_redirect(request, pk, group=None, user=None):
    # Handler when redirecting to a group view.
    # Permissions to check if the user may view the group at all is at middleware level.
    return group_redirect_handler(request=request, group=group, user=user, autoadd=False)


def group_redirect_handler(request, group, user, autoadd=None):
    # Redirects to a group.
    try:
        target = auth.get_group_url(group)

        if group.public and user.is_authenticated() and autoadd:
            # Only public groups may be automatically joined.
            auth.groupsub_get_or_create(user=user, usergroup=group)

        return redirect(target)

    except Exception as exc:
        messages.error(request, "Group error: %s" % exc)
        return redirect(reverse("home"))

@auth.group_access
def group_info(request, pk, group=None, user=None):
    template_name = "group_info.html"

    site = models.Site.objects.get(id=settings.SITE_ID)

    print (site.domain)
    # Current group permissions
    perms = GroupPerm.objects.filter(usergroup=group).select_related("user")

    context = dict(perms=perms, target=group)
    return render(request, template_name, context)


def group_list(request):
    """
    Generates the list of groups.
    """
    template_name = "group_list.html"

    user, group = request.user, request.group

    if user.is_authenticated():
        # See if the user has permission to the current group.
        curr_perm = GroupPerm.objects.filter(user=user, usergroup=group, role=GroupPerm.ADMIN).first()

        # Find all group subscriptions for the user.
        sub_list = GroupSub.objects.filter(user=user)

        # Turn into a dictionary for fast lookup
        sub_map = dict((s.usergroup, s) for s in sub_list)

        # Add all groups to the filtering condition
        ids = [sub.usergroup_id for sub in sub_list]
        cond = Q(public=True) | Q(id__in=ids)
    else:
        curr_perm = None
        sub_map = {}
        cond = Q(public=True)

    # Decorate the current group with a role attribute.
    group.curr_perm = curr_perm

    groups = UserGroup.objects.filter(cond).select_related("owner")

    paginator = query.ExtendedPaginator(request,
                                        object_list=groups, per_page=25)
    page = paginator.curr_page()

    for g in page.object_list:
        g.editable = (g.owner == user)
        sub = sub_map.get(g)
        g.subscription = sub.get_type_display() if sub else "None"

    context = dict(page=page, public=page.object_list, group=group)

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

    # Add the recent votes
    recent_votes = query.recent_votes()

    context = dict(page=page, posts=page.object_list,
                   recent_votes=recent_votes, q=q)

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
            PostView.objects.create(ip=ip, post=post, date=auth.now())
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
    write_access_check = auth.thread_write_access(user=user, root=post)

    def decorator(p):
        # Each post needs to carry information on its status relative to the user.
        p.editable = write_access_check(user=user, post=p)
        p.has_vote = p.id in post.upvotes
        p.has_bookmark = p.id in post.bookmarks
        return p

    # Decorate all posts in the thread.
    thread = list(map(decorator, thread))

    # Decorate the main post as well.
    post = decorator(post)

    # Store answers in a separate list for simpler access.
    post.answers = filter(lambda p: p.type == models.Post.ANSWER, thread)
    post.answers = list(post.answers)

    # Comments will be stored in a dictionary for fast access.
    comment_list = filter(lambda pc: pc.type == models.Post.COMMENT, thread)

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

