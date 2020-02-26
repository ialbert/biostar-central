import datetime
import logging
import hashlib
import urllib.parse
from django.template import loader
from django.conf import settings
from django.contrib import messages
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import F, Q
from django.utils.timezone import utc
from django.core.cache import cache
from biostar.accounts.models import Profile, Logger
from . import util
from .const import *
from .models import Post, Vote, PostView, Subscription

User = get_user_model()

logger = logging.getLogger("engine")


def get_votes(user, root):
    store = {
        Vote.BOOKMARK: set(),
        Vote.UP: set()
    }

    # Collect all the votes for the user.
    if user.is_authenticated:
        votes = Vote.objects.filter(post__root=root, author=user).values_list("type", "post__id")

        for vote_type, post_id, in votes:
            store.setdefault(vote_type, set()).add(post_id)

    return store


def gravatar_url(email, style='mp', size=80):
    hash_num = hashlib.md5(email).hexdigest()

    url = "https://secure.gravatar.com/avatar/%s?" % hash_num
    url += urllib.parse.urlencode({
        's': str(size),
        'd': style,
    }
    )
    return url


def get_users_str():
    """
    Return comma separated string of username used for autocomplete
    """

    cache_days = 5
    cache_secs = 60 * 60 * 24 * cache_days

    users_str = cache.get(USERS_CACHE_KEY)
    if users_str is None:
        users_str = ','.join(User.objects.all().values_list('username', flat=True))
        cache.set(USERS_CACHE_KEY, users_str, cache_secs)

    return users_str



def gravatar(user, size=80):
    if not user or user.is_anonymous:
        email = 'anon@biostars.org'.encode('utf8')
        return gravatar_url(email=email)

    email = user.email if user.is_authenticated else ''
    email = email.encode('utf8')

    if user.is_anonymous or not user.profile.is_valid:
        # Removes spammy images for suspended users
        email = 'suspended@biostars.org'.encode('utf8')

        style = settings.GRAVATAR_ICON or "monsterid"
    elif user.profile.is_moderator:
        style = settings.GRAVATAR_ICON or "robohash"
    elif user.profile.score > 100:
        style = settings.GRAVATAR_ICON or "retro"
    elif user.profile.score > 0:
        style = settings.GRAVATAR_ICON or "identicon"
    else:
        style = settings.GRAVATAR_ICON or "mp"

    return gravatar_url(email=email, style=style, size=size)


def walk_down_thread(parent, collect=set()):
    """
    Recursively walk down a thread of posts starting from target
    """

    # Stop condition: post does not have a root or parent.
    if (parent is None) or (parent.parent is None) or (parent.root is None):
        return collect

    # Get all children for this post, excluding itself.
    children = Post.objects.filter(parent=parent).exclude(uid=parent.uid)

    for child in children:
        # Add child to list
        collect.add(child)
        # Get all children belonging to the current child.
        walk_down_thread(parent=child, collect=collect)

    return collect


def create_subscription(post, user, sub_type=None, update=False):
    """
    Creates subscription to a post. Returns a list of subscriptions.
    """
    subs = Subscription.objects.filter(post=post.root, user=user)
    sub = subs.first()

    default = Subscription.TYPE_MAP.get(user.profile.message_prefs,
                                        Subscription.LOCAL_MESSAGE)
    sub_type = sub_type or (sub.type if sub else None) or default

    # Ensure the sub type is not set to something wrote
    if sub and update:
        # Update an existing subscription
        sub.type = sub_type
        sub.save()
    else:
        # Drop all existing subscriptions for the user by default.
        subs.delete()
        Subscription.objects.create(post=post.root, user=user, type=sub_type)

    # Recompute post subscription.
    subs_count = Subscription.objects.filter(post=post.root).exclude(type=Profile.NO_MESSAGES).count()

    # Update root subscription counts.
    Post.objects.filter(pk=post.root.pk).update(subs_count=subs_count)


def is_suspended(user):
    if user.is_authenticated and user.profile.state in (Profile.BANNED, Profile.SUSPENDED, Profile.SPAMMER):
        return True

    return False


def post_tree(user, root):
    """
    Populates a tree that contains all posts in the thread.

    Answers sorted before comments.
    """

    # Get all posts that belong to post root.
    query = Post.objects.filter(root=root).exclude(pk=root.id)

    query = query.select_related("lastedit_user__profile", "author__profile", "root__author__profile")

    is_moderator = user.is_authenticated and user.profile.is_moderator

    # Only moderators
    if not is_moderator:
        query = query.exclude(status=Post.DELETED)
        # query = query.exclude(spam=Post.SPAM)

    # Apply the sort order to all posts in thread.
    thread = query.order_by("type", "-accept_count", "-vote_count", "creation_date")

    # Gather votes by the current user.
    votes = get_votes(user=user, root=root)

    # Shortcuts to each storage.
    bookmarks, upvotes = votes[Vote.BOOKMARK], votes[Vote.UP]

    # Build comments tree.
    comment_tree = dict()

    def decorate(post):
        # Mutates the elements! Not worth creating copies.
        if post.is_comment:
            comment_tree.setdefault(post.parent_id, []).append(post)
        post.has_bookmark = int(post.id in bookmarks)
        post.has_upvote = int(post.id in upvotes)
        post.can_accept = not post.is_toplevel and (user == post.root.author or (user.is_authenticated and user.profile.is_moderator))
        post.can_moderate = user.is_authenticated and user.profile.is_moderator
        post.is_editable = user.is_authenticated and (user == post.author or (user.is_authenticated and user.profile.is_moderator))
        return post

    # Decorate the objects for easier access
    thread = list(map(decorate, thread))

    # Decorate the root post
    root = decorate(root)

    # Select the answers from the thread.
    answers = [p for p in thread if p.type == Post.ANSWER]

    return root, comment_tree, answers, thread


def update_post_views(post, request, minutes=settings.POST_VIEW_MINUTES):
    "Views are updated per user session"

    # Extract the IP number from the request.
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    # 'localhost' is not a valid ip address.
    ip1 = '' if ip1.lower() == 'localhost' else ip1
    ip2 = '' if ip2.lower() == 'localhost' else ip2
    ip = ip1 or ip2 or '0.0.0.0'

    now = util.now()
    since = now - datetime.timedelta(minutes=minutes)

    # One view per time interval from each IP address.
    if not PostView.objects.filter(ip=ip, post=post, date__gt=since).exists():
        # Update the last time
        PostView.objects.create(ip=ip, post=post, date=now)
        Post.objects.filter(pk=post.pk).update(view_count=F('view_count') + 1)
    return post


@transaction.atomic
def apply_vote(post, user, vote_type):
    vote = Vote.objects.filter(author=user, post=post, type=vote_type).first()

    if vote:
        msg = "%s removed" % vote.get_type_display()
        change = -1
        vote.delete()
    else:
        change = +1
        vote = Vote.objects.create(author=user, post=post, type=vote_type)
        msg = "%s added" % vote.get_type_display()

    if post.author == user:
        # Author making the change
        change = 0
        return msg, vote, change

    # Fetch user score
    score = len(Vote.objects.filter(post__author=post.author).exclude(author=post.author))
    # Update the user score.
    Profile.objects.filter(user=post.author).update(score=score)

    # Calculate counts for the current post
    votes = list(Vote.objects.filter(post=post).exclude(author=post.author))
    vote_count = len(votes)
    bookcount = len(list(filter(lambda v: v.type == Vote.BOOKMARK, votes)))
    accept_count = len(list(filter(lambda v: v.type == Vote.ACCEPT, votes)))

    # Increment the post vote count.
    Post.objects.filter(uid=post.uid).update(vote_count=vote_count)

    # The thread vote count represents all votes in a thread
    Post.objects.filter(uid=post.root.uid).update(thread_votecount=F('thread_votecount') + change)

    # Increment the bookmark count.
    if vote_type == Vote.BOOKMARK:
        Post.objects.filter(uid=post.uid).update(book_count=bookcount)

    # Handle accepted vote.
    if vote_type == Vote.ACCEPT:
        Post.objects.filter(uid=post.uid).update(accept_count=accept_count)
        Post.objects.filter(uid=post.root.uid).update(accept_count=F('accept_count') + change)

    return msg, vote, change


def log_action(user=None, action=Logger.MODERATING, log_text=''):
    # Create a logger object
    Logger.objects.create(user=user, action=action, log_text=log_text)
    return


def delete_post(post, request):
    # Delete marks a post deleted but does not remove it.

    # Posts with children or older than some value can only be deleted not removed
    # The children of a post.
    children = Post.objects.filter(parent_id=post.id).exclude(pk=post.id)

    # The condition where post can only be deleted.
    delete_only = children or post.age_in_days > 7 or post.vote_count > 1 or (post.author != request.user)

    if delete_only:
        # Deleted posts can be undeleted by re-opening them.
        Post.objects.filter(uid=post.uid).update(status=Post.DELETED)
        url = post.root.get_absolute_url()
        messages.success(request, "Deleted post: %s" % post.title)
    # Remove post from the database with no trace.
    else:
        # This will remove the post. Redirect depends on the level of the post.
        url = "/" if post.is_toplevel else post.root.get_absolute_url()
        post.delete()
        messages.success(request, "Removed post: %s" % post.title)

    # Recompute answers count
    if post.type == Post.ANSWER:
        answer_count = Post.objects.filter(root=post.root, type=Post.ANSWER).count()
        Post.objects.filter(pk=post.parent_id).update(answer_count=answer_count)

    reply_count = Post.objects.filter(root=post.root).count()

    Post.objects.filter(pk=post.root.id).update(reply_count=reply_count)
    log_action(user=request.user, log_text=f"Deleted post={post.uid}")

    return url


def handle_spam_post(post, user):
    url = post.get_absolute_url()

    # # Ban new users that post spam.
    # if post.author.profile.low_rep:
    #     post.author.profile.state = Profile.BANNED

    post.author.profile.state = Profile.SPAMMER
    post.author.profile.save()

    # Label all posts by this users as spam.
    Post.objects.filter(author=post.author).update(spam=Post.SPAM, status=Post.OFFTOPIC)
    log_action(user=user, log_text=f"Reported post={post.uid} as spam.")
    return url


def moderate_post(request, action, post, offtopic='', comment=None, dupes=[], pid=None):
    root = post.root
    user = request.user
    now = datetime.datetime.utcnow().replace(tzinfo=utc)
    url = post.get_absolute_url()

    if action == BUMP_POST:
        Post.objects.filter(uid=post.uid).update(lastedit_date=now, rank=now.timestamp())
        messages.success(request, "Post bumped")
        log_action(user=user, log_text=f"Bumped post={post.uid}")
        return url

    if action == OPEN_POST:
        Post.objects.filter(uid=post.uid).update(status=Post.OPEN, spam=Post.NOT_SPAM)
        messages.success(request, f"Opened post: {post.title}")
        log_action(user=user, log_text=f"Opened post={post.uid}")
        return url

    if action == DELETE:
        return delete_post(post=post, request=request)

    if action == MOVE_ANSWER:
        Post.objects.filter(uid=post.uid).update(type=Post.ANSWER)
        log_action(user=user, log_text=f"Moved post={post.uid} to answer. ")
        return url

    if action == REPORT_SPAM:
        return handle_spam_post(post=post, user=user)

    if pid:
        parent = Post.objects.filter(uid=pid).first() or post.root
        Post.objects.filter(uid=post.uid).update(type=Post.COMMENT, parent=parent)
        Post.objects.filter(uid=root.uid).update(reply_count=F("answer_count") - 1)
        messages.success(request, "Moved answer to comment")
        log_action(user=user, log_text=f"Moved post={post.uid} to comment.")
        return url

    if dupes and len(''.join(dupes)):
        # Load comment explaining post off topic label.
        tmpl = loader.get_template("messages/duplicate_posts.md")
        context = dict(user=post.author, dupes=dupes, comment=comment)
        content = tmpl.render(context)

        post = Post.objects.create(content=content, type=Post.ANSWER, parent=post, root=post.root, author=user)
        url = post.get_absolute_url()
        messages.success(request, "Marked duplicated post as off topic.")
        log_action(user=user, log_text=f"Marked post={post.uid} as duplicate.")
        return url

    if comment:
        # Load comment explaining post closure.
        tmpl = loader.get_template("messages/off_topic.md")
        context = dict(user=post.author, comment=comment)
        content = tmpl.render(context)

        # Load answer explaining post being off topic.
        post = Post.objects.create(content=content, type=Post.ANSWER, parent=post, root=post.root, author=user)
        url = post.get_absolute_url()
        messages.success(request, "Marked the post as off topic.")
        log_action(user=user, log_text=f"Marked post={post.uid} as off topic.")
        return url

    return url
