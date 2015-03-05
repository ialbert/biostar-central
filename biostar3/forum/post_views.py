from __future__ import absolute_import, division, print_function, unicode_literals

# Python modules.
from collections import OrderedDict, defaultdict

# Django specific modules.
from django.shortcuts import render, redirect
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse

# Biostar specific local modules.
from . import models, query, search, auth
from .models import Vote, Post, PostView

# Get custom user model.
User = get_user_model()


def post_list(request):
    template_name = "post_list.html"
    posts = query.get_toplevel_posts(user=request.user, group=request.group)
    page = query.get_page(request, posts, per_page=settings.POSTS_PER_PAGE)

    # Add the recent votes
    recent_votes = query.recent_votes()
    html_title = "Post List"
    context = dict(page=page, posts=page.object_list, recent_votes=recent_votes, html_title=html_title)

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

    page = query.get_page(request, posts, per_page=settings.POSTS_PER_PAGE)

    # Add the recent votes
    recent_votes = query.recent_votes()
    html_title = "Post List"
    context = dict(page=page, posts=page.object_list, recent_votes=recent_votes, html_title=html_title)

    return render(request, template_name, context)


def post_view(request, pk):
    """
    Generates the page that contains a full thread.
    """
    user = request.user
    template_name = "post_detail.html"

    # Tries to get the post.
    post = Post.objects.filter(pk=pk).first()

    if not post:
        # Post does not exist.
        messages.error(request, "This post does not exist. Perhaps it has been deleted.")
        return redirect("home")

    if not auth.read_access_post(user=user, post=post):
        # Post exists but may not be read by the user.
        messages.error(request, "This post my not be accessed by this user.")
        return redirect("home")

    if not post.is_toplevel:
        # Post is not at top level. Redirect and and scroll the page to the right anchor.
        return redirect(post.get_absolute_url())

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
    thread = map(decorator, thread)

    # Decorate the main post as well.
    post = decorator(post)

    # Store answers in a separate list for simpler access.
    post.answers = filter(lambda p: p.type == models.Post.ANSWER, thread)

    # Comments will be stored in a dictionary for fast access.
    comment_list = filter(lambda pc: pc.type == models.Post.COMMENT, thread)

    # Collect comments into a dictionary keyed by the parent id with
    # comments as a value list
    post.comments = OrderedDict()
    for comment in comment_list:
        post.comments.setdefault(comment.parent.id, []).append(comment)


    if user.is_authenticated():
        # This is for testing only. Keeps adding comments to posts on the page.
        import random, faker

        f = faker.Factory.create()
        u = random.choice(User.objects.all())
        parent = random.choice(thread + [post])
        text = f.bs()
        comment = models.Post.objects.create(type=models.Post.COMMENT, parent=parent, content=text, author=u)

    # Add object to the context.
    html_title = post.title
    context = dict(post=post, html_title=html_title)

    return render(request, template_name, context)

