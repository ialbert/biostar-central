import logging
import os
from datetime import timedelta
from functools import wraps

from django.conf import settings
from django.contrib import messages
from django.contrib.auth import get_user_model
from django.contrib.auth.decorators import login_required
from django.core.cache import cache
from django.core.paginator import Paginator
from django.db.models import Count, Q
from django.http import Http404
from django.shortcuts import render, redirect, reverse
from django.views.decorators.csrf import ensure_csrf_cookie
from taggit.models import Tag
from biostar.planet.models import Blog, BlogPost
from biostar.accounts.models import Profile
from biostar.forum import forms, auth, tasks, util, search, models, moderate
from biostar.forum.const import *

from biostar.forum.models import Post, Vote, Badge, Subscription, Log
from biostar.utils.decorators import is_moderator, check_params, reset_count, is_staff, authenticated

User = get_user_model()

logger = logging.getLogger('engine')

RATELIMIT_KEY = settings.RATELIMIT_KEY

CREATE_PARAMS = {'title', 'tag_val'}
CREATE_PARAMS.update(ALLOWED_PARAMS)

# Valid post values as they correspond to database post types.
POST_TYPE = dict(
    question=Post.QUESTION,
    jobs=Post.JOB,
    tutorials=Post.TUTORIAL,
    forum=Post.FORUM,
    planet=Post.BLOG,
    tools=Post.TOOL,
    news=Post.NEWS,
    pages=Post.PAGE,
    herald=Post.HERALD
)

LIMIT_MAP = dict(
    all=0,
    today=1,
    week=7,
    month=30,
    year=365
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


class CachedPaginator(Paginator):
    """
    Paginator that caches the count call.
    """

    # Time to live for the cache, in seconds
    TTL = 300

    def __init__(self, cache_key='', ttl=None, *args, **kwargs):
        self.cache_key = cache_key

        # May not contain spaces
        self.cache_key = ''.join(self.cache_key.split())

        self.ttl = self.TTL

        super(CachedPaginator, self).__init__(*args, **kwargs)

    @property
    def count(self):

        if self.cache_key:
            # See if it is access the cache
            value = cache.get(self.cache_key)
            if value is None:
                value = super(CachedPaginator, self).count
                # logger.debug(f'setting the cache for "{self.cache_key}"')
                cache.set(self.cache_key, value, self.ttl)
        else:
            value = super(CachedPaginator, self).count

        return value


def apply_sort(posts, limit=None, order=None):
    # Apply post ordering.
    if ORDER_MAPPER.get(order):
        ordering = ORDER_MAPPER.get(order)
        posts = posts.order_by(ordering)
    else:
        posts = posts.order_by('-rank')

    days = LIMIT_MAP.get(limit, 0)
    # Apply time limit if required.
    if days:
        delta = util.now() - timedelta(days=days)
        posts = posts.filter(lastedit_date__gt=delta)

    # Select related information used during rendering.
    posts = posts.select_related("root").select_related("author__profile", "lastedit_user__profile")

    return posts


@reset_count(key='spam_count')
def get_spam(request):
    posts = Post.objects.filter(spam=Post.SPAM)
    return posts


def get_posts(request, topic=""):
    """
    Generates a post list on a topic.
    """

    user = request.user
    # Topics are case insensitive.
    topic = topic or LATEST
    topic = topic.lower()

    # Detect known post types.
    post_type = POST_TYPE.get(topic)

    # Get all open top level posts.
    posts = Post.objects.filter(is_toplevel=True, status=Post.OPEN)

    # Pass the above query when asking for latest
    if topic == LATEST:
        pass
    # Filter for various post types.
    elif post_type is not None:
        posts = posts.filter(type=post_type)

    elif topic == SHOW_SPAM and (user.is_authenticated and user.profile.is_moderator):
        posts = get_spam(request)

    elif topic == OPEN:
        posts = posts.filter(type=Post.QUESTION, answer_count=0)

    elif topic == BOOKMARKS and user.is_authenticated:
        posts = Post.objects.filter(votes__author=user, votes__type=Vote.BOOKMARK)

    elif topic == FOLLOWING and user.is_authenticated:
        posts = posts.filter(subs__user=user).exclude(subs__type=Subscription.NO_MESSAGES)

    elif topic == MYPOSTS and user.is_authenticated:
        # Show users all of their posts ( deleted, spam, or quarantined )
        posts = Post.objects.filter(author=user)

    elif topic == MYVOTES and user.is_authenticated:
        posts = posts.filter(votes__post__author=user)

    elif topic == MYTAGS and user.is_authenticated:
        tags = map(lambda t: t.lower(), user.profile.my_tags.split(","))
        posts = posts.filter(tags__name__in=tags).distinct()

    else:
        posts = Post.objects.none()

    return posts


@check_params(allowed=ALLOWED_PARAMS)
def post_search(request):
    query = request.GET.get('query', '')
    length = len(query.replace(" ", ""))
    page = int(request.GET.get('page', 1))
    order = request.GET.get('order', 'relevance')

    mapper = dict(relevance=None, date=['lastedit_date'])
    sortedby = mapper.get(order)

    if length < settings.SEARCH_CHAR_MIN:
        messages.error(request, "Enter more characters before preforming search.")
        return redirect(reverse('post_list'))

    # Reverse sort when ordering by date.
    revsort = order == 'date'
    results, indexed = search.perform_search(query=query, page=page, reverse=revsort, sortedby=sortedby)

    context = dict(results=results, query=query, indexed=indexed, order=order)

    return render(request, "search/search_results.html", context=context)


@check_params(allowed=ALLOWED_PARAMS)
def pages(request, fname):
    # Add markdown file extension to markdown
    infile = f"{fname}.md"
    # Look for this file in static root.
    doc = os.path.join(settings.STATIC_ROOT, "forum", infile)

    if not os.path.exists(doc):
        messages.error(request, "File does not exist.")
        return redirect("post_list")

    admins = User.objects.filter(profile__role=Profile.MANAGER)
    mods = User.objects.filter(profile__role=Profile.MODERATOR).exclude(id__in=admins)
    admins = admins.prefetch_related("profile").order_by("-profile__score")
    mods = mods.prefetch_related("profile").order_by("-profile__score")
    context = dict(file_path=doc, tab=fname, admins=admins, mods=mods)

    return render(request, 'pages.html', context=context)


@is_moderator
def mark_spam(request, uid):
    """
    Mark post as spam.
    """

    # Trigger post
    post = Post.objects.filter(uid=uid).first()

    # A restore parameter sent toggles spam off.
    state = False if request.GET.get("restore") else True

    # Apply the toggle.
    if post:
        moderate.toggle_spam(request, post)
    else:
        messages.error(request, "Post does not seem to exist")

    # Was spam actually
    if post.is_spam:
        return redirect(reverse('post_topic', kwargs=dict(topic='spam')))
    else:
        return redirect('/')


@is_moderator
def release_quar(request, uid):
    """
    Release quarantined post to the public.
    """
    post = Post.objects.filter(uid=uid).first()
    if not post:
        messages.error(request, "Post does not exist.")
        return redirect('/')

    # Bump the score by one is the user does not get quarantined again.
    # Tells the system user has gained antibodies!
    if post.author.profile.low_rep:
        post.author.profile.bump_over_threshold()

    Post.objects.filter(uid=uid).update(spam=Post.NOT_SPAM)

    return redirect('/')


def post_list(request, topic=None, tag="", cutoff=None, ordering=None):
    """
    Post listing. Filters, orders and paginates posts based on GET parameters.
    """

    # Parse the GET parameters for filtering information
    page = request.GET.get('page', 1)
    order = request.GET.get("order", ordering) or 'rank'
    topic = topic or request.GET.get("type", LATEST) or LATEST
    limit = request.GET.get("limit", "all") or "all"

    if tag:
        # Get all open top level posts.
        posts = Post.objects.filter(is_toplevel=True, status=Post.OPEN, tags__name__iexact=tag)
        cache_key = ''
    else:
        # Get posts available to users.
        posts = get_posts(request=request, topic=topic)
        # Create the cache key only with latest topic
        cache_key = f"{LATEST}-{order}-{limit}" if topic is LATEST else ''

    posts = apply_sort(posts, limit=limit, order=order)

    # Institute a cutoff
    if cutoff:
        posts = posts[:cutoff]

    # Filter for any empty strings
    paginator = CachedPaginator(cache_key=cache_key, object_list=posts, per_page=settings.POSTS_PER_PAGE)

    # Apply the post paging.
    posts = paginator.get_page(page)

    return posts


@check_params(allowed=ALLOWED_PARAMS)
@ensure_csrf_cookie
def latest(request):
    """
    Show latest post listing.
    """

    posts = post_list(request, topic=LATEST)

    context = dict(posts=posts, tab=LATEST)

    return render(request, template_name="post_list.html", context=context)


@check_params(allowed=ALLOWED_PARAMS)
@ensure_csrf_cookie
def post_tags(request, tag):
    """
    Show list of posts belonging to one post.
    """
    posts = post_list(request, tag=tag, cutoff=settings.POSTS_PER_PAGE)
    # Clear tags if no posts are found for it
    tag = tag if posts else ''
    context = dict(posts=posts, tag=tag)

    return render(request, template_name="post_tags.html", context=context)


@check_params(allowed=ALLOWED_PARAMS)
@ensure_csrf_cookie
def post_topic(request, topic):
    """
    Show list of posts of a given type
    """

    # Set the cache key based on order and limit
    posts = post_list(request, topic=topic, cutoff=1000)

    # Clear topic if there are no posts.
    topic = topic if posts else ''
    context = dict(posts=posts, topic=topic, tab=topic)

    return render(request, template_name="post_topic.html", context=context)


@ensure_csrf_cookie
@authenticated
def bookmarks(request):
    """
    Show posts bookmarked by user.
    """

    posts = post_list(request, topic=BOOKMARKS, ordering=VOTE_DATE)

    # Order by vote date.
    context = dict(posts=posts, topic=BOOKMARKS, tab=BOOKMARKS)
    return render(request, template_name="user_bookmarks.html", context=context)


@ensure_csrf_cookie
@authenticated
def mytags(request):
    posts = post_list(request, topic=MYTAGS)

    context = dict(posts=posts, topic=MYTAGS, tab=MYTAGS)
    return render(request, template_name="user_mytags.html", context=context)


@ensure_csrf_cookie
@authenticated
def myposts(request):
    """
    Show posts by user
    """
    posts = post_list(request, topic=MYPOSTS)

    context = dict(posts=posts, topic=MYPOSTS, tab=MYPOSTS)

    return render(request, template_name="user_myposts.html", context=context)


@ensure_csrf_cookie
@authenticated
def following(request):
    """
    Show posts followed by user.
    """
    posts = post_list(request, topic=FOLLOWING)

    context = dict(posts=posts, topic=FOLLOWING, tab=FOLLOWING)
    return render(request, template_name="user_following.html", context=context)


@authenticated
@reset_count(key="vote_count")
def myvotes(request):
    """
    Show posts by user that received votes
    """
    page = request.GET.get('page', 1)

    votes = Vote.objects.filter(post__author=request.user).select_related('post', 'post__root',
                                                                          'author__profile').order_by("-date")
    # Create the paginator
    paginator = CachedPaginator(object_list=votes,
                                per_page=settings.POSTS_PER_PAGE)

    # Apply the votes paging.
    votes = paginator.get_page(page)

    context = dict(votes=votes, page=page, tab='myvotes')
    return render(request, template_name="user_votes.html", context=context)


@check_params(allowed=ALLOWED_PARAMS)
def tags_list(request):
    """
    Show posts by user
    """
    page = request.GET.get('page', 1)
    query = request.GET.get('query', '')

    count = Count('post', filter=Q(post__is_toplevel=True))

    db_query = Q(name__icontains=query) if query else Q()
    cache_key = '' if query else TAGS_CACHE_KEY

    tags = Tag.objects.annotate(nitems=count).filter(db_query)
    tags = tags.order_by('-nitems')

    # Create the paginator
    paginator = CachedPaginator(cache_key=cache_key,
                                object_list=tags,
                                per_page=settings.POSTS_PER_PAGE)

    # Apply the votes paging.
    tags = paginator.get_page(page)

    context = dict(tags=tags, tab='tags', query=query)

    return render(request, 'tags_list.html', context=context)


@check_params(allowed=ALLOWED_PARAMS)
def community_list(request):
    page = request.GET.get("page", 1)
    ordering = request.GET.get("order", "visit")
    limit_to = request.GET.get("limit", "time")
    query = request.GET.get('query', '')
    query = query.replace("'", "").replace('"', '').strip()
    days = LIMIT_MAP.get(limit_to, 0)

    users = User.objects.select_related("profile")

    if days:
        delta = util.now() - timedelta(days=days)
        users = users.filter(profile__last_login__gt=delta)

    if query and len(query) > 2:
        db_query = Q(profile__name__icontains=query) | Q(profile__uid__icontains=query) | \
                   Q(username__icontains=query) | Q(email__icontains=query)
        users = users.filter(db_query)

    order = ORDER_MAPPER.get(ordering, "visit")
    users = users.filter(profile__state__in=[Profile.NEW, Profile.TRUSTED])
    users = users.order_by(order)

    # Create the paginator (six users per row)
    paginator = CachedPaginator(object_list=users, per_page=60)
    users = paginator.get_page(page)
    context = dict(tab="community", users=users, query=query, order=ordering, limit=limit_to)

    return render(request, "community_list.html", context=context)


@check_params(allowed=ALLOWED_PARAMS)
def badge_list(request):
    badges = Badge.objects.annotate(count=Count("award")).order_by('-count')
    context = dict(badges=badges)
    return render(request, "badge_list.html", context=context)


@check_params(allowed=ALLOWED_PARAMS)
def badge_view(request, uid):
    badge = Badge.objects.filter(uid=uid).annotate(count=Count("award")).first()
    target = request.GET.get('user')
    page = request.GET.get('page', 1)

    if not badge:
        messages.error(request, f"Badge with id={uid} does not exist.")
        return redirect(reverse("badge_list"))

    awards = badge.award_set.all().order_by("-date")

    if target:
        user = User.objects.filter(profile__uid=target).first()
        awards = awards.filter(user=user)

    awards = awards.prefetch_related("user", "user__profile", "post", "post__root")
    paginator = Paginator(object_list=awards, per_page=settings.POSTS_PER_PAGE)

    awards = paginator.get_page(page)
    context = dict(awards=awards, badge=badge)

    return render(request, "badge_view.html", context=context)


@check_params(allowed=ALLOWED_PARAMS)
@ensure_csrf_cookie
def post_view(request, uid):
    "Return a detailed view for specific post"

    # Get the post.
    post = Post.objects.filter(uid=uid).select_related('root').first()
    user = request.user
    if not post:
        messages.error(request, "Post does not exist.")
        return redirect("post_list")

    # Redirect to post view
    if not post.is_toplevel:
        return redirect(post.get_absolute_url())

    # Return 404 when post is spam and user is anonymous.
    if post.is_spam and user.is_anonymous:
        raise Http404("Post does not exist.")

    # Form used for answers
    form = forms.PostShortForm(user=request.user, post=post)

    if request.method == "POST":
        form = forms.PostShortForm(data=request.POST, ptype=Post.ANSWER, user=request.user, post=post)
        if form.is_valid():
            answer = form.save()
            return redirect(answer.get_absolute_url())

        messages.error(request, form.errors)

    # Build the comment tree .
    root, comment_tree, answers, thread = auth.post_tree(user=request.user, root=post.root)

    # Bump post views.
    models.update_post_views(post=post, request=request, timeout=settings.POST_VIEW_TIMEOUT)

    context = dict(post=root, tree=comment_tree, form=form, answers=answers)

    return render(request, "post_view.html", context=context)


@check_params(allowed=CREATE_PARAMS)
@login_required
def new_post(request):
    """
    Creates a new post
    """
    title = request.GET.get('title', '')
    tag_val = request.GET.get('tag_val', '')
    tag_val = ','.join(tag_val.split())
    initial = dict(title=title, tag_val=tag_val)
    content = ''
    author = request.user
    form = forms.PostLongForm(user=request.user, initial=initial)
    if request.method == "POST":

        form = forms.PostLongForm(data=request.POST, user=request.user, initial=initial)
        tag_val = form.data.get('tag_val')
        content = form.data.get('content', '')
        if form.is_valid():
            # Create a new post by user
            title = form.cleaned_data.get('title')
            content = form.cleaned_data.get("content")
            ptype = form.cleaned_data.get('post_type')
            tag_val = form.cleaned_data.get('tag_val')
            post = auth.create_post(title=title, content=content, ptype=ptype, tag_val=tag_val, author=author,
                                    request=request)

            tasks.created_post.spool(pid=post.id)

            return redirect(post.get_absolute_url())

    # Action url for the form is the current view
    action_url = reverse("post_create")
    context = dict(form=form, tab="new", tag_val=tag_val, action_url=action_url,
                   content=content)

    return render(request, "new_post.html", context=context)


@reset_count(key='mod_count')
@check_params(allowed=ALLOWED_PARAMS)
@login_required
def view_logs(request):
    LIMIT = 100

    if request.user.profile.is_moderator:
        logs = Log.objects.all().order_by("-id")[:LIMIT]

    else:
        logs = Log.objects.filter(pk=0)

    logs = logs.select_related("user", "post", "post__root", "user__profile", "target", "target__profile",
                               "post__author")

    context = dict(logs=logs)

    return render(request, "view_logs.html", context=context)


@is_staff
def merge_profile(request):
    """
    Merge two profiles into one.
    """

    user = request.user
    form = forms.MergeProfiles(user=user)

    if request.method == 'POST':
        form = forms.MergeProfiles(user=user, data=request.POST)

        if form.is_valid():
            merged = form.save()
            messages.success(request, "Merged profiles")
            return redirect(reverse('user_profile', kwargs=dict(uid=merged.profile.uid)))

    context = dict(form=form)
    return render(request, "accounts/merge_profile.html", context=context)


def error(request):
    """
    Checking error propagation and logging
    """
    1 / 0
