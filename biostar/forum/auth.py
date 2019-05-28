import datetime
import logging

from django.conf import settings
from django.contrib import messages
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import F
from django.utils.timezone import utc

from biostar.accounts.models import Profile
from . import util
from .const import *
from .models import Post, Vote, Subscription, PostView

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
    query = query.select_related("lastedit_user__profile", "root__lastedit_user__profile",
                                 "root__author__profile", "author__profile")

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
        post.is_editable = user.is_authenticated and (user == post.author or user.profile.is_moderator)
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
    else:
        # Update the various counts only if the user is different.
        Profile.objects.filter(user=post.author).update(score=F('score') + change)

        # Increment the post vote count.
        Post.objects.filter(uid=post.uid).update(vote_count=F('vote_count') + change)

        # The thread vote count represents all votes in a thread
        Post.objects.filter(uid=post.root.uid).update(thread_votecount=F('thread_votecount') + change)

        # Increment the bookmark count.
        if vote_type == Vote.BOOKMARK:
            Post.objects.filter(uid=post.uid).update(book_count=F('book_count') + change)

        # Handle accepted vote.
        if vote_type == Vote.ACCEPT:
            Post.objects.filter(uid=post.uid).update(accept_count=F('accept_count') + change)
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


def moderate_post(request, action, post, comment=None, dupes=[], pid=None):
    root = post.root
    user = request.user
    now = datetime.datetime.utcnow().replace(tzinfo=utc)
    url = post.root.get_absolute_url()

    if action == BUMP_POST:
        Post.objects.filter(uid=post.uid).update(lastedit_date=now, lastedit_user=request.user)
        messages.success(request, "Post bumped")
        return url

    if action == MOD_OPEN:
        Post.objects.filter(uid=post.uid).update(status=Post.OPEN)
        messages.success(request, f"Opened post: {post.title}")
        return url

    if action == DELETE:
        return delete_post(post=post, request=request)

    if action == TOGGLE_ACCEPT:
        # Recompute accept count for post
        change = -1 if post.has_accepted else + 1
        Post.objects.filter(uid=post.uid).update(accept_count=F("accept_count") + change)
        Post.objects.filter(uid=root.uid).update(accept_count=F("accept_count") + change)

        return url

    if action == MOVE_TO_ANSWER:
        Post.objects.filter(uid=post.uid).update(type=Post.ANSWER, parent=post.root)
        Post.objects.filter(uid=root.uid).update(reply_count=F("reply_count") + 1)
        messages.success(request, "Moved comment to answer")
        return url

    if action == MOVE_TO_COMMENT or pid:
        parent = Post.objects.filter(uid=pid).first() or post.root
        Post.objects.filter(uid=post.uid).update(type=Post.COMMENT, parent=parent)
        Post.objects.filter(uid=root.uid).update(reply_count=F("reply_count") - 1)
        messages.success(request, "Moved answer to comment")
        return url

    if dupes:
        Post.objects.filter(uid=post.uid).update(status=Post.CLOSED)
        html = util.render(name="default_messages/duplicate_posts.html", user=post.author, dupes=dupes,
                           comment=comment, posts=Post.objects.filter(uid=post.uid))
        content = util.strip_tags(html)
        # Create a comment to the post
        modpost = Post.objects.create(content=content, type=Post.COMMENT, parent=post, author=user)
        Post.objects.filter(uid=modpost.uid).update(html=html)
        return url

    messages.error(request, "Invalid moderation action given")
    return url


def create_post(author, content, post_type, title="Title", tag_val="tag1, tag2", parent=None, root=None,):
    "Used to create posts across apps"

    post = Post.objects.create(title=title, content=content, tag_val=tag_val,
                               author=author, type=post_type, parent=parent, root=root)

    # Trigger notifications for subscribers and mentioned users
    # async or synchronously

    return post
