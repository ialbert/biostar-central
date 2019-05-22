import logging
from datetime import timedelta

from django.conf import settings
from django.contrib import messages
from django.contrib.auth import get_user_model
from django.contrib.auth.decorators import login_required
from django.core.paginator import Paginator
from django.db.models import Count
from django.shortcuts import render, redirect

from biostar.forum import forms, auth, tasks, util
from biostar.forum.const import *
from biostar.forum.models import Post, Vote, Subscription, Badge
from biostar.utils.decorators import ajax_error, ajax_error_wrapper, ajax_success, object_exists
from biostar.utils.shortcuts import reverse

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
    views="-view_count",
    replies="-reply_count",
    votes="-thread_votecount",
    visit='-profile__last_login',
    reputation='-profile__score',
    joined='-profile__date_joined',
    activity='-profile__date_joined'

)


def get_posts(user, show="latest", tag="", order="rank", limit=None):
    """
    Generates a post list on a topic.
    """
    # Topics are case insensitive.
    topic = show.lower()

    # Detect known post types.
    post_type = POST_TYPE_MAPPER.get(topic)
    # Determines how to start the query.
    if post_type:
        query = Post.objects.filter(type=post_type)
    elif topic == OPEN:
        query = Post.objects.filter(type=Post.QUESTION, reply_count=0)
    elif topic == BOOKMARKS and user.is_authenticated:
        query = Post.objects.filter(votes__author=user, votes__type=Vote.BOOKMARK)
    elif topic == FOLLOWING and user.is_authenticated:
        query = Post.objects.exclude(subs__type=Subscription.NO_MESSAGES).filter(subs__user=user)
    elif topic == MYPOSTS and user.is_authenticated:
        query = Post.objects.filter(author=user)
    elif topic == MYVOTES and user.is_authenticated:
        # TODO: switching to votes
        # votes_query = Vote.objects.filter(post__author=user).exclude(author=user)
        # query = votes_query.values("post")
        query = Post.objects.filter(votes__post__author=user).exclude(votes__author=user)
    else:
        query = Post.objects.filter(type__in=Post.TOP_LEVEL)

    # Filter by tags if specified.
    if tag:
        query = query.filter(tag_val__iregex=tag)

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

    # Filter deleted items for non subscribed users
    cond = user.is_authenticated and user.profile.is_moderator
    query = query if cond else query.exclude(status=Post.DELETED)
    # Select related information used during rendering.
    query = query.prefetch_related("root", "author__profile", "lastedit_user__profile", "thread_users__profile")

    return query


def post_list(request, show=None):
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

    # Create the paginator
    paginator = Paginator(posts, settings.POSTS_PER_PAGE)

    # Apply the post paging.
    posts = paginator.get_page(page)

    # Set the active tab.
    tab = show or tag or "latest"

    # Fill in context.
    context = dict(posts=posts, tab=tab, tag=tag, order=order, limit=limit)

    # Render the page.
    return render(request, template_name="post_list.html", context=context)


def authenticated(func):
    def _wrapper_(request, **kwargs):
        if request.user.is_anonymous:
            messages.error(request, "You need to be logged in to view this page.")
        return func(request, **kwargs)

    return _wrapper_


def latest(request):
    show = request.GET.get("type", "") or LATEST
    return post_list(request, show=show)


@authenticated
def myvotes(request):
    """
    Show posts by user that received votes
    """
    return post_list(request, show=MYVOTES)


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


def community_list(request):
    users = User.objects.select_related("profile")
    page = request.GET.get("page", 1)
    ordering = request.GET.get("order", "visit")
    limit_to = request.GET.get("limit", "time")
    days = LIMIT_MAP.get(limit_to, 0)

    if days:
        delta = util.now() - timedelta(days=days)
        users = users.filter(profile__last_login__gt=delta)

    order = ORDER_MAPPER.get(ordering, "visit")
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

@ajax_error_wrapper(method="POST")
def comment_box(request):
    """
    Returns the comment box
    """


@ajax_error_wrapper(method="POST")
def ajax_vote(request):
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

    msg, vote = auth.apply_vote(post=post, user=user, vote_type=vote_type)

    jquery.poshytip.js


@object_exists(klass=Post)
def post_view(request, uid):
    "Return a detailed view for specific post"

    # Form used for answers
    form = forms.PostShortForm()

    # Get the post.
    post = Post.objects.filter(uid=uid).first()

    # Establish the root for the post.
    root = post if post.is_toplevel else post.root

    # auth.update_post_views(post=obj, request=request)

    # Build the comment tree .
    root, comment_tree, answers, thread = auth.post_tree(user=request.user, root=root)

    context = dict(post=post, tree=comment_tree, form=form, answers=answers)

    return render(request, "post_view.html", context=context)


@object_exists(klass=Post)
def post_answer(request, uid):
    """
    Process an answer with form
    """
    # Get the post.
    root = Post.objects.filter(uid=uid).first()
    redir = reverse("post_view", request=request, kwargs=dict(uid=root.uid))
    if request.method == "POST":
        form = forms.PostShortForm(data=request.POST)
        if form.is_valid():
            answer = form.save(author=request.user)
            if tasks.HAS_UWSGI:
                tasks.created_post(pid=answer.id)

            # Anchor location to recently created answer
            redir = redir + "#" + answer.uid

    return redirect(redir)


def ajax_test(request):
    """
    Creates a commment on a top level post.
    """
    msg="OK"
    print (f"HeRe= {request.POST} ")
    return ajax_success(msg=msg)


@login_required
def create_comment(request, uid):
    user = request.user
    post = Post.objects.filter(uid=uid).first()

    if request.method == "POST":
        form = forms.CommentForm(data=request.POST)
        if form.is_valid():
            content = form.cleaned_data['content']
            comment = auth.create_post(parent=post.parent, author=user, content=content, post_type=Post.COMMENT)
            context = dict(comment=comment)
            return render(request, "widgets/render_comment.html", context=context)

        messages.error(request, f"Error adding comment:{form.errors}")
    else:
        initial = dict(parent_uid=post.uid, content="")
        form = forms.CommentForm(initial=initial)

    context = dict(post=post, form=form)
    return render(request, "widgets/create_comment.html", context=context)


@object_exists(klass=Post)
@login_required
def subs_action(request, uid):
    # Post actions are being taken on
    post = Post.objects.filter(uid=uid).first()
    user = request.user

    if request.method == "POST" and user.is_authenticated:
        form = forms.SubsForm(data=request.POST, post=post, user=user)

        if form.is_valid():
            sub = form.save()
            msg = f"Updated Subscription to : {sub.get_type_display()}"
            messages.success(request, msg)

            if tasks.HAS_UWSGI:
                tasks.added_sub(sid=sub.id)

    return redirect(reverse("post_view", kwargs=dict(uid=post.uid)))


@login_required
def post_create(request, template="post_create.html", url="post_view",
                extra_context={}, filter_func=lambda x: x):
    "Make a new post"

    # Filter function ( filter_func ) is used to filter choices from the form
    # between sites.
    form = forms.PostLongForm(filter_func=filter_func)

    if request.method == "POST":
        form = forms.PostLongForm(data=request.POST, filter_func=filter_func)
        if form.is_valid():
            # Create a new post by user
            post = form.save(author=request.user)
            if tasks.HAS_UWSGI:
                tasks.created_post(pid=post.id)
            return redirect(reverse(url, request=request, kwargs=dict(uid=post.uid)))
    context = dict(form=form, tab="new", action_url=reverse("post_create"))
    context.update(extra_context)

    return render(request, template, context=context)


@object_exists(klass=Post)
@login_required
def post_moderate(request, uid):
    user = request.user
    post = Post.objects.filter(uid=uid).first()

    if request.method == "POST":

        form = forms.PostModForm(post=post, data=request.POST, user=user, request=request)

        if form.is_valid():

            action = form.cleaned_data["action"]
            duplicate = form.cleaned_data["dupe"]
            pid = form.cleaned_data.get("pid", "")
            redir = auth.moderate_post(post=post, request=request, action=action, dupes=duplicate, pid=pid)
            if tasks.HAS_UWSGI:
                tasks.moderated_post(pid=post.id)
            return redirect(redir)
        else:
            messages.error(request, "Invalid moderation error.")
            return redirect(reverse("post_view", kwargs=dict(uid=uid)))
    else:
        form = forms.PostModForm(post=post, user=user, request=request)

    context = dict(form=form, post=post)
    return render(request, "post_moderate.html", context)


@object_exists(klass=Post)
@login_required
def edit_post(request, uid):
    "Edit an existing post"

    post = Post.objects.filter(uid=uid).first()
    if post.is_toplevel:
        template, edit_form = "post_create.html", forms.PostLongForm
    else:
        template, edit_form = "shortpost_edit.html", forms.PostShortForm

    user = request.user
    form = edit_form(post=post, user=user)
    if request.method == "POST":
        form = edit_form(post=post, data=request.POST, user=user)
        if form.is_valid():
            form.save(edit=True)
            messages.success(request, f"Edited :{post.title}")
            location = reverse("post_view", kwargs=dict(uid=post.root.uid)) + "#" + post.uid
            return redirect(location)

    context = dict(form=form, post=post, action_url=reverse("post_edit", kwargs=dict(uid=uid)),
                   extra_tab="active", extra_tab_name="Edit Post")

    return render(request, template, context)
