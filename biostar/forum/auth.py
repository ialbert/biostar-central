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

from biostar.accounts.models import Profile
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
    hash = hashlib.md5(email).hexdigest()

    url = "https://secure.gravatar.com/avatar/%s?" % hash
    url += urllib.parse.urlencode({
        's': str(size),
        'd': style,
    }
    )
    return url


def gravatar(user, size=80):

    if not user or user.is_anonymous:
        email = 'anon@biostars.org'.encode('utf8')
        return gravatar_url(email=email)

    email = user.email if user.is_authenticated else ''
    email = email.encode('utf8')

    if user.is_anonymous or user.profile.is_suspended or user.profile.is_banned:
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


def walk_down_thread(parent, collect=[], is_root=True):
    """
    Recursively walk up a thread of posts starting from target
    """

    # Stop condition 1: post does not have a root or parent.
    if (parent is None) or (parent.parent is None) or (parent.root is None):
        return collect

    # Get all children for this post
    if is_root:
        children = Post.objects.filter(root=parent).exclude(uid=parent.uid)
    else:
        children = Post.objects.filter(parent=parent).exclude(uid=parent.uid)

    for child in children:
        # Add child to list
        collect.append(child)
        # Get all children belonging to the current child.
        walk_down_thread(parent=child, collect=collect, is_root=is_root)

    return collect


def create_subscription(post, user, sub_type=None, delete_exisiting=True):
    """
    Creates subscription to a post. Returns a list of subscriptions.
    """
    # Drop all existing subscriptions for the user by default.
    if delete_exisiting:
        Subscription.objects.filter(post=post.root, user=user).delete()
        # Create new subscription to the user.
        Subscription.objects.create(post=post.root, user=user, type=sub_type)
    # Update an existing subscription type.
    else:
        sub, created = Subscription.objects.get_or_create(post=post.root, user=user)
        Subscription.objects.filter(pk=sub.pk).update(type=sub_type)

    # Recompute post subscription.
    subs_count = Subscription.objects.filter(post=post.root).exclude(type=Profile.NO_MESSAGES).count()

    # Update root subscription counts.
    Post.objects.filter(pk=post.root.pk).update(subs_count=subs_count)


def is_suspended(user):
    return user.is_authenticated and user.profile.state in (Profile.BANNED, Profile.SUSPENDED)


def post_tree(user, root):
    """
    Populates a tree that contains all posts in the thread.

    Answers sorted before comments.
    """

    # Get all posts that belong to post root.
    query = Post.objects.filter(root=root).exclude(pk=root.id)

    # Add all related objects.
    query = query.select_related("lastedit_user__profile", "lastedit_user", "root__lastedit_user",
                                 "root__lastedit_user__profile", "root__author__profile", "author__profile")

    is_moderator = user.is_authenticated and user.profile.is_moderator

    # Only moderators
    if not is_moderator:
        query = query.exclude(status=Post.DELETED)

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
    answers = [ p for p in thread if p.type==Post.ANSWER ]

    return root, comment_tree,  answers, thread


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

    return url


def moderate_post(request, action, post, offtopic='', comment=None, dupes=[], pid=None):
    root = post.root
    user = request.user
    now = datetime.datetime.utcnow().replace(tzinfo=utc)
    url = post.get_absolute_url()

    if action == BUMP_POST:
        Post.objects.filter(uid=post.uid).update(lastedit_date=now, rank=now.timestamp(), last_contributor=request.user)
        messages.success(request, "Post bumped")
        return url

    if action == OPEN_POST:
        Post.objects.filter(uid=post.uid).update(status=Post.OPEN)
        messages.success(request, f"Opened post: {post.title}")
        return url

    if action == DELETE:
        return delete_post(post=post, request=request)

    if action == MOVE_ANSWER:
        Post.objects.filter(uid=post.uid).update(type=Post.ANSWER)
        return url

    if action == LOCK:
        Post.objects.filter(uid=post.uid).update(status=Post.LOCKED)
        messages.success(request, f"Locked post: {post.title}")
        return url

    if action == CLOSE:
        Post.objects.filter(uid=post.uid).update(status=Post.CLOSED)
        messages.success(request, f"Closed post: {post.title}")
        return url

    if pid:
        parent = Post.objects.filter(uid=pid).first() or post.root
        Post.objects.filter(uid=post.uid).update(type=Post.COMMENT, parent=parent)
        Post.objects.filter(uid=root.uid).update(reply_count=F("answer_count") - 1)
        messages.success(request, "Moved answer to comment")
        return url

    if offtopic:
        # Load comment explaining post closure.
        tmpl = loader.get_template("messages/off_topic.md")
        context = dict(user=post.author, comment=offtopic)
        content = tmpl.render(context)

        Post.objects.filter(uid=post.uid).update(status=Post.OFFTOPIC)
        # Load answer explaining post being off topic.
        post = Post.objects.create(content=content, type=Post.ANSWER, parent=post, author=user)
        url = post.get_absolute_url()
        messages.success(request, "Marked the post as off topic.")

        return url

    if dupes:
        # Load comment explaining post off topic label.
        tmpl = loader.get_template("messages/duplicate_posts.md")
        context = dict(user=post.author, dupes=dupes, comment=comment)
        content = tmpl.render(context)

        Post.objects.filter(uid=post.uid).update(status=Post.OFFTOPIC)
        post = Post.objects.create(content=content, type=Post.COMMENT, parent=post, author=user)
        url = post.get_absolute_url()
        messages.success(request, "Closed duplicated post.")

        return url

    return url