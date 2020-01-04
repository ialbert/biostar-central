import logging
from datetime import timedelta
from functools import wraps

from whoosh.searching import Results
from django.conf import settings
from django.contrib import messages
from django.contrib.auth import get_user_model
from django.contrib.auth.decorators import login_required
from django.views.decorators.csrf import ensure_csrf_cookie
from django.core.paginator import Paginator
from django.db.models import Count, Q
from taggit.models import Tag
from django.shortcuts import render, redirect, reverse

from biostar.accounts.models import Profile
from . import forms, auth, tasks, util, search
from .const import *
from .models import Post, Vote, Badge


User = get_user_model()

logger = logging.getLogger('engine')

# Valid post values as they correspond to database post types.
POST_TYPE_MAPPER = dict(
    question=Post.QUESTION,
    job=Post.JOB,
    tutorial=Post.TUTORIAL,
    forum=Post.FORUM,
    blog=Post.BLOG,
    tool=Post.TOOL,
    news=Post.NEWS
)

LIMIT_MAP = dict(
    all=0,
    today=1,
    week=7,
    month=30,
    year=365
)
# Valid order values value as they correspond to database ordering fields.
ORDER_MAPPER = dict(
    rank="-rank",
    tagged='-tagged',
    views="-view_count",
    replies="-reply_count",
    votes="-thread_votecount",
    visit='-profile__last_login',
    reputation='-profile__score',
    joined='-profile__date_joined',
    activity='-profile__date_joined'


)


def post_exists(func):
    """
    Ensure uid passed to view function exists.
    """
    @wraps(func)
    def _wrapper_(request, **kwargs):
        uid = kwargs.get('uid')
        post = Post.objects.filter(uid=uid).exists()
        if not post:
            messages.error(request, "Post does not exist.")
            return redirect(reverse("post_list"))
        return func(request, **kwargs)
    return _wrapper_


def policy(request):
    context = dict()
    return render(request, template_name="widgets/policy.html", context=context)


def get_posts(user, show="latest", tag="", order="rank", limit=None):
    """
    Generates a post list on a topic.
    """
    # Topics are case insensitive.
    topic = show.lower()

    # Detect known post types.
    post_type = POST_TYPE_MAPPER.get(topic)
    # Determines how to start the preform_search.
    if post_type:
        query = Post.objects.filter(type=post_type)
    elif topic == OPEN:
        query = Post.objects.filter(type=Post.QUESTION, answer_count=0)
    elif topic == BOOKMARKS and user.is_authenticated:
        query = Post.objects.filter(votes__author=user, votes__type=Vote.BOOKMARK)
    elif topic == FOLLOWING and user.is_authenticated:
        query = Post.objects.filter(subs__user=user)
    elif topic == MYPOSTS and user.is_authenticated:
        query = Post.objects.filter(author=user)
    elif topic == MYVOTES and user.is_authenticated:
        query = Post.objects.filter(votes__post__author=user)
    elif topic == MYTAGS and user.is_authenticated:
        tags = user.profile.my_tags.split(",")
        query = Post.objects.filter(tags__name__in=tags)
    else:
        query = Post.objects.filter(is_toplevel=True)

    # Filter by tags if specified.
    if tag:
        query = query.filter(tags__name=tag.lower())

    # Apply post ordering.
    if ORDER_MAPPER.get(order):
        ordering = ORDER_MAPPER.get(order)
        query = query.order_by(ordering)
    else:
        query = query.order_by("-rank")

    days = LIMIT_MAP.get(limit, 0)
    # Apply time limit if required.
    if days:
        delta = util.now() - timedelta(days=days)
        query = query.filter(lastedit_date__gt=delta)

    # Filter deleted items for anonymous and non-moderators.
    if user.is_anonymous or (user.is_authenticated and not user.profile.is_moderator):
        query = query.exclude(status=Post.DELETED)

    # Select related information used during rendering.
    query = query.prefetch_related("root", "author__profile", "lastedit_user__profile", "thread_users__profile")

    return query


def post_search(request):

    query = request.GET.get('query', '')
    page = int(request.GET.get('page', 1))

    if not query:

        return redirect(reverse('post_list'))

    # Preform search on indexed posts.
    results = search.preform_whoosh_search(query=query, page=page, per_page=settings.SEARCH_RESULTS_PER_PAGE)

    question_flag = Post.QUESTION
    context = dict(results=results, query=query, question_flag=question_flag, stop_words=','.join(search.STOP))

    return render(request, template_name="widgets/post_results.html", context=context)


@ensure_csrf_cookie
def post_list(request, show=None, extra_context=dict()):
    """
    Post listing. Filters, orders and paginates posts based on GET parameters.
    """
    # The user performing the request.
    user = request.user

    # Parse the GET parameters for filtering information
    page = request.GET.get('page', 1)
    tag = request.GET.get("tag", "")
    order = request.GET.get("order", "rank")
    show = show or request.GET.get("type", "")
    limit = request.GET.get("limit", "all")

    # Get posts available to users.
    posts = get_posts(user=user, show=show, tag=tag, order=order, limit=limit)
    posts = posts.exclude(Q(root=None) | Q(parent=None))
    # Create the paginator
    paginator = Paginator(posts, settings.POSTS_PER_PAGE)

    # Apply the post paging.
    posts = paginator.get_page(page)

    # Set the active tab.
    tab = tag or show or "latest"

    # Fill in context.
    context = dict(posts=posts, tab=tab, tag=tag, order=order, type=show, limit=limit)
    context.update(extra_context)
    # Render the page.
    return render(request, template_name="post_list.html", context=context)


def latest(request):
    show = request.GET.get("type", "") or LATEST
    return post_list(request, show=show)


def authenticated(func):
    def _wrapper_(request, **kwargs):
        if request.user.is_anonymous:
            messages.error(request, "You need to be logged in to view this page.")
        return func(request, **kwargs)
    return _wrapper_


@authenticated
def myvotes(request):
    """
    Show posts by user that received votes
    """
    page = request.GET.get('page', 1)
    votes = Vote.objects.filter(post__author=request.user).order_by("-date")
    # Create the paginator
    paginator = Paginator(votes, settings.POSTS_PER_PAGE)

    # Apply the votes paging.
    votes = paginator.get_page(votes)

    context = dict(votes=votes, page=page, tab='myvotes')
    return render(request, template_name="votes_list.html", context=context)


def tags_list(request):
    """
    Show posts by user
    """
    page = request.GET.get('page', 1)

    count = Count('post', filter=Q(post__type__in=Post.TOP_LEVEL))
    tags = Tag.objects.annotate(nitems=count)
    tags = tags.order_by('-nitems')
    # Create the paginator
    paginator = Paginator(tags, 150)

    # Apply the votes paging.
    tags = paginator.get_page(page)

    context = dict(tags=tags, tab='tags')

    return render(request, 'tags_list.html', context=context)


@authenticated
def myposts(request):
    """
    Show posts by user
    """
    return post_list(request, show=MYPOSTS)


@authenticated
def following(request):
    """
    Show posts followed by user
    """
    return post_list(request, show=FOLLOWING)


@authenticated
def bookmarks(request):
    """
    Show posts bookmarked by user
    """
    return post_list(request, show=BOOKMARKS)


@authenticated
def mytags(request):

    return post_list(request=request, show=MYTAGS)


def community_list(request):
    users = User.objects.select_related("profile")
    page = request.GET.get("page", 1)
    ordering = request.GET.get("order", "visit")
    limit_to = request.GET.get("limit", "time")
    query = request.GET.get('query', '')
    days = LIMIT_MAP.get(limit_to, 0)

    if days:
        delta = util.now() - timedelta(days=days)
        users = users.filter(profile__last_login__gt=delta)

    if query and len(query) > 2:
        db_query = Q(email__in=query) | Q(profile__name__contains=query) | Q(profile__uid__contains=query) | \
                   Q(username__contains=query) | Q(profile__name__in=query) | Q(email=query) |Q(email__contains=query) |\
                   Q(profile__uid__contains=query)
        users = users.filter(db_query)

    order = ORDER_MAPPER.get(ordering, "visit")
    users = users.exclude(profile__state__in=[Profile.SUSPENDED, Profile.BANNED])
    users = users.order_by(order)

    paginator = Paginator(users, settings.USERS_PER_PAGE)
    users = paginator.get_page(page)
    context = dict(tab="community", users=users, order=ordering, limit=limit_to)

    return render(request, "community_list.html", context=context)


def badge_list(request):
    badges = Badge.objects.annotate(count=Count("award"))
    context = dict(badges=badges)
    return render(request, "badge_list.html", context=context)


def badge_view(request, uid):
    badge = Badge.objects.filter(uid=uid).annotate(count=Count("award")).first()

    if not badge:
        messages.error(request, f"Badge with id={uid} does not exist.")
        return redirect(reverse("badge_list"))

    awards = badge.award_set.order_by("-pk")[:100]
    awards = awards.prefetch_related("user", "user__profile", "post", "post__root")
    context = dict(awards=awards, badge=badge)

    return render(request, "badge_view.html", context=context)


# def tags_list(request):

#    context = dict(extra_tab="active", extra_tab_name="All Tags")
#    return render(request, "tags_list.html", context=context)


@ensure_csrf_cookie
@post_exists
def post_view(request, uid):
    "Return a detailed view for specific post"

    # Get the post.
    post = Post.objects.filter(uid=uid).first()

    auth.update_post_views(post=post, request=request)
    if not post.is_toplevel:
        return redirect(post.get_absolute_url())

    # Form used for answers
    form = forms.PostShortForm(user=request.user, post=post)

    if request.method == "POST":

        form = forms.PostShortForm(data=request.POST, user=request.user, post=post)
        if form.is_valid():
            author = request.user
            content = form.cleaned_data.get("content")
            # Create answer to root
            answer = Post.objects.create(title=post.title, parent=post, author=author,
                                         content=content, type=Post.ANSWER, root=post.root)
            tasks.created_post.spool(pid=answer.id)
            return redirect(answer.get_absolute_url())
        messages.error(request, form.errors)

    # Build the comment tree .
    root, comment_tree, answers, thread = auth.post_tree(user=request.user, root=post.root)

    context = dict(post=root, tree=comment_tree, form=form, answers=answers)

    return render(request, "post_view.html", context=context)


def new_comment(request, uid):
    user = request.user
    post = Post.objects.filter(uid=uid).first()

    if request.method == "POST":
        form = forms.PostShortForm(data=request.POST, recaptcha=False, user=user)
        if form.is_valid():
            content = form.cleaned_data['content']
            comment = Post.objects.create(parent=post, author=user, content=content, type=Post.COMMENT)
            return redirect(comment.get_absolute_url())
        messages.error(request, f"Error adding comment:{form.errors}")
    else:
        initial = dict(parent_uid=post.uid, content="")
        form = forms.PostShortForm(initial=initial, recaptcha=False, user=user)

    context = dict(post=post, form=form, user=user)

    return render(request, "new_comment.html", context=context)


@login_required
def new_post(request):
    """
    Creates a new post
    """
    form = forms.PostLongForm(user=request.user)
    author = request.user

    if request.method == "POST":

        form = forms.PostLongForm(data=request.POST, user=request.user)
        if form.is_valid():
            # Create a new post by user
            title = form.cleaned_data.get('title')
            content = form.cleaned_data.get("content")
            post_type = form.cleaned_data.get('post_type')
            tag_val = form.cleaned_data.get('tag_val')
            post = Post.objects.create(title=title, content=content, type=post_type,
                                       tag_val=tag_val, author=author)

            tasks.created_post.spool(pid=post.id)

            return redirect(post.get_absolute_url())

    # Action url for the form is the current view
    action_url = reverse("post_create")
    content = request.POST.get('content', '')
    context = dict(form=form, tab="new", action_url=action_url, content=content)

    return render(request, "new_post.html", context=context)


@post_exists
@login_required
def post_moderate(request, uid):
    """Used to make dispaly post moderate form given a post request."""

    user = request.user
    post = Post.objects.filter(uid=uid).first()

    if not user.profile.is_moderator:
        messages.error(request, 'You need to be a moderator to perform this action.')
        return redirect(post.get_absolute_url())

    if post.status == Post.LOCKED:
        messages.error(request, 'Only admins can moderate locked posts.')
        return redirect(post.get_absolute_url())

    if request.method == "POST":
        form = forms.PostModForm(post=post, data=request.POST, user=user, request=request)

        if form.is_valid():
            action = form.cleaned_data.get('action')
            dupe = form.cleaned_data.get('dupe', '').split(",")
            dupe_comment = form.cleaned_data.get('comment')
            mod_uid = form.cleaned_data.get('pid')
            offtopic = form.cleaned_data.get('offtopic')
            redir = auth.moderate_post(post=post, request=request, action=action, comment=dupe_comment,
                                       dupes=dupe, pid=mod_uid, offtopic=offtopic)
            return redirect(redir)
        else:
            messages.error(request, "Invalid action")
            return redirect(reverse("post_view", kwargs=dict(uid=uid)))
    else:
        form = forms.PostModForm(post=post, user=user, request=request)

    context = dict(form=form, post=post)
    return render(request, "post_moderate.html", context)


@post_exists
@login_required
def edit_post(request, uid):
    """
    Edit an existing post"
    """
    post = Post.objects.filter(uid=uid).first()
    action_url = reverse("post_edit", kwargs=dict(uid=post.uid))
    user = request.user
    initial = dict(content=post.content, title=post.title, tag_val=post.tag_val, post_type=post.type)

    form = forms.PostLongForm(post=post, initial=initial, user=user)

    if request.method == "POST":
        form = forms.PostLongForm(post=post, initial=initial, data=request.POST, user=user)
        if form.is_valid():
            form.edit()
            messages.success(request, f"Edited :{post.title}")
            return redirect(post.get_absolute_url())

    context = dict(form=form, post=post, action_url=action_url, form_title="Edit post", content=post.content)

    return render(request, "new_post.html", context)
