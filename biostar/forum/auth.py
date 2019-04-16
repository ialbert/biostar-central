
import bleach
import datetime
import logging
import re
import mistune
from itertools import chain

from django.contrib import messages
from django.utils.timezone import utc
from django.db.models import F, Q
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction

from biostar.message import tasks
from biostar.forum.awards import ALL_AWARDS
from biostar.accounts.models import Profile
from biostar.utils.shortcuts import reverse
from .models import Post, Vote, Subscription, PostView, Award, Badge
from . import util
from .const import *

User = get_user_model()


logger = logging.getLogger("engine")


def get_votes(user, thread):

    store = {Vote.BOOKMARK: set(), Vote.UP:set()}

    if user.is_authenticated:
        votes = Vote.objects.filter(post__in=thread, author=user).values_list("post__id", "type")

        for post_id, vote_type in votes:
            store.setdefault(vote_type, set()).add(post_id)

    return store


def my_posts(target, request):

    user = request.user
    if user.is_anonymous or target.is_anonymous:
        return Post.objects.filter(author=target).exclude(status=Post.DELETED)

    query = Post.objects.filter(author=target)
    query = query if (user.profile.is_moderator or user == target) else query.exclude(status=Post.DELETED)

    return query


def build_obj_tree(request, obj):

    # Populate the object to build a tree that contains all posts in the thread.
    # Answers sorted before comments.
    user = request.user
    query = Post.objects.filter(root=obj)
    query = query if user.is_authenticated and user.profile.is_moderator else query.exclude(status=Post.DELETED)
    thread = query.order_by("type", "-has_accepted", "-vote_count", "creation_date")

    thread = thread.select_related("lastedit_user__profile", "root__author__profile",
                                   "author__profile")
    # Gather votes
    votes = get_votes(user=user, thread=thread)

    # Shortcuts to each storage.
    bookmarks = votes[Vote.BOOKMARK]
    upvotes = votes[Vote.UP]
    # Build comments tree.
    comment_tree = dict()

    def decorate(posts):
        # Can the current user accept answers
        # TODO: use annotate.

        for post in posts:
            if post.is_comment:
                comment_tree.setdefault(post.parent_id, []).append(post)

            post.has_bookmark = post.id in bookmarks
            post.has_upvote = post.id in upvotes
            post.is_editable = user.is_authenticated and (user == post.author or user.profile.is_moderator)

    answers = thread.filter(type=Post.ANSWER)

    # Decorate the objects for easier access
    decorate(chain(thread, answers))

    return comment_tree, answers, thread


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


def create_sub(post,  user, sub_type=None):
    "Creates a subscription of a user to a post"

    root = post.root
    sub = Subscription.objects.filter(post=root, user=user).first()
    date = datetime.datetime.utcnow().replace(tzinfo=utc)

    # Update an existing sub with new type.
    if sub:
        Subscription.objects.filter(pk=sub.pk).update(type=sub_type)
        # The sub is being changed to "No message"
        if sub_type == Subscription.NO_MESSAGES:
            Post.objects.filter(pk=root.pk).update(subs_count=F('subs_count') - 1)

    # Create a new sub object
    else:
        sub = Subscription.objects.create(post=root, user=user, type=sub_type, date=date)
        # Increase the subscription count of the root.
        if sub_type != Subscription.NO_MESSAGES:
            Post.objects.filter(pk=root.pk).update(subs_count=F('subs_count') + 1)

    return sub


def trigger_vote(vote_type, post, change):
    Post.objects.get_all(uid=post.uid).update(vote_count=F('vote_count') + change)

    if vote_type == Vote.BOOKMARK:
        # Apply the vote
        Post.objects.get_all(uid=post.uid).update(book_count=F('book_count') + change)

    elif vote_type == Vote.ACCEPT:

        if change > 0:
            # There does not seem to be a negation operator for F objects.
            Post.objects.get_all(uid=post.uid).update(has_accepted=True)
            Post.objects.get_all(uid=post.root.uid).update(has_accepted=True)
        else:
            Post.objects.get_all(uid=post.uid).update(has_accepted=False)
            accepted_siblings = Post.objects.get_all(root=post.root, has_accepted=True).exclude(pk=post.root_id).count()

            # Only set root as not accepted if there are no accepted siblings
            if accepted_siblings == 0:
                Post.objects.get_all(uid=post.root.uid).update(has_accepted=False)
    else:
        thread_query = Post.objects.filter(status=Post.OPEN, root=post.root)

        reply_count = thread_query.exclude(uid=post.parent.uid).filter(type=Post.ANSWER).count()
        thread_score = thread_query.exclude(uid=post.root.uid).count()
        Post.objects.get_all(root=post.root).update(thread_votecount=F('thread_votecount') + change)
        Post.objects.filter(parent=post.parent).update(reply_count=reply_count)
        Post.objects.filter(root=post.root).update(thread_score=thread_score)


@transaction.atomic
def preform_vote(post, user, vote_type, uid=None):

    vote = Vote.objects.filter(author=user, post=post, type=vote_type).first()

    if vote:
        msg = "%s removed" % vote.get_type_display()
        change = -1
        vote.delete()
    else:
        change = +1
        vote = Vote.objects.create(author=user, post=post, type=vote_type, uid=uid)
        msg = "%s added" % vote.get_type_display()

    if post.author != user:
        # Update the user reputation only if the author is different.
        Profile.objects.filter(user=post.author).update(score=F('score') + change)

    # The thread vote count represents all votes in a thread
    Post.objects.get_all(uid=post.root.uid).update(thread_votecount=F('thread_votecount') + change)

    trigger_vote(vote_type=vote_type, post=post, change=change)

    return msg, vote


def create_post_from_json(json_dict):

    root_uid = json_dict.get("root_id")
    parent_uid = json_dict.get("parent_id", None)

    lastedit_user_uid = json_dict.get("lastedit_user_id")
    author_uid = json_dict.get("author_id")

    author = User.objects.filter(profile__uid=author_uid).first()
    lastedit_user = User.objects.filter(profile__uid=lastedit_user_uid).first()

    root = Post.objects.filter(uid=root_uid).first()
    parent = Post.objects.filter(uid=parent_uid).first()

    creation_date = json_dict.get("creation_date")
    lastedit_date = json_dict.get("lastedit_date")

    title = json_dict.get("title")
    has_accepted = json_dict.get("has_accepted", False)
    type = json_dict.get("type")
    status = json_dict.get("status", Post.OPEN)
    content = util.strip_tags(json_dict.get("text", ""))
    html = json_dict.get("html", "")
    tag_val = json_dict.get("tag_val")

    #reply_count = json_dict.get("reply_count", 0)
    #thread_score = json_dict.get("thread_score", 0)
    #vote_count = json_dict.get("vote_count", 0)
    view_count = json_dict.get("view_count", 0)

    uid = json_dict.get("id")
    post = Post.objects.filter(uid=uid)
    if post.exists() or status == Post.DELETED:
        logger.error(f"Post with uid={uid} already exists or status is deleted.")
        return post.first()

    post = Post.objects.create(uid=uid, author=author, lastedit_user=lastedit_user,
                               root=root, parent=parent, creation_date=creation_date,
                               lastedit_date=lastedit_date, title=title, has_accepted=has_accepted,
                               type=type, status=status, content=content, html=html, tag_val=tag_val,
                               view_count=view_count)
    # Trigger another save
    #post.add_tags(post.tag_val)

    logger.info(f"Created post.uid={post.uid}")

    return post


def parse_mentioned_users(content):

    # Any word preceded by a @ is considered a user handler.
    handler_pattern = "\@[^\s]+"
    # Drop leading @
    users_list = set(x[1:] for x in re.findall(handler_pattern, content))

    return User.objects.filter(username__in=users_list)


def parse_html(text):
    "Sanitize text and expand links to match content"

    # This will collect the objects that could be embedded
    mentioned_users = parse_mentioned_users(content=text)

    html = mistune.markdown(text)

    # embed the objects
    for user in mentioned_users:
        url = reverse("user_profile", kwargs=dict(uid=user.profile.uid))
        handler = f"@{user.username}"
        emb_patt = f'<a href="{url}">{handler}</a>'
        html = html.replace(handler, emb_patt)

    return html


def delete_post(post, request):
    # Delete marks a post deleted but does not remove it.
    # Remove means to delete the post from the database with no trace.

    # Posts with children or older than some value can only be deleted not removed
    # The children of a post.
    children = Post.objects.filter(parent_id=post.id).exclude(pk=post.id)

    # The condition where post can only be deleted.
    delete_only = children or post.age_in_days > 7 or post.vote_count > 1 or (post.author != request.user)

    if delete_only:
        # Deleted posts can be undeleted by re-opening them.
        Post.objects.get_all(uid=post.uid).update(status=Post.DELETED)
        url = post.root.get_absolute_url()
        messages.success(request, "Deleted post: %s" % post.title)
    else:
        # This will remove the post. Redirect depends on the level of the post.
        url = "/" if post.is_toplevel else post.root.get_absolute_url()
        post.delete()
        messages.success(request, "Removed post: %s" % post.title)

    # Recompute post reply count
    if post.type == Post.ANSWER:
        reply_count = Post.objects.filter(parent=post.parent, type=Post.ANSWER, status=Post.OPEN).count()
        Post.objects.filter(pk=post.parent_id).update(reply_count=reply_count)

    thread_score = Post.objects.filter(type=Post.ANSWER, root=post.root, status=Post.OPEN).count()
    Post.objects.filter(pk=post.root_id).update(thread_score=thread_score)

    return url


def moderate_post(request, action, post, comment=None, dupes=[]):

    root = post.root
    user = request.user
    now = datetime.datetime.utcnow().replace(tzinfo=utc)
    url = post.root.get_absolute_url()

    if action == BUMP_POST:
        Post.objects.get_all(uid=post.uid).update(lastedit_date=now, lastedit_user=request.user)
        messages.success(request, "Post bumped")
        return url

    if action == MOD_OPEN:
        Post.objects.get_all(uid=post.uid).update(status=Post.OPEN)
        messages.success(request, f"Opened post: {post.title}")
        return url

    if action == DELETE:
        return delete_post(post=post, request=request)

    if action == CROSSPOST:
        content = util.render(name="messages/crossposted.html", user=post.author, comment=comment, posts=post)
        # Create a comment to the post
        Post.objects.create(content=content, type=Post.COMMENT, html=content, parent=post, author=user)
        return url

    if action == TOGGLE_ACCEPT:
        root_has_accepted = Post.objects.get_all(root=root, type=Post.ANSWER, has_accepted=True).count()
        Post.objects.get_all(uid=post.uid).update(has_accepted=not post.has_accepted)
        Post.objects.get_all(uid=root.uid).update(has_accepted=root_has_accepted)
        return url

    if action == MOVE_TO_ANSWER:
        Post.objects.get_all(uid=post.uid).update(type=Post.ANSWER, parent=post.root, reply_count=F("reply_count") + 1)
        Post.objects.get_all(uid=root.uid).update(reply_count=F("reply_count") + 1)
        messages.success(request, "Moved comment to answer")
        return url

    if action == MOVE_TO_COMMENT:
        Post.objects.get_all(uid=post.uid).update(type=Post.COMMENT, parent=post.root, reply_count=F("reply_count") - 1)
        Post.objects.get_all(uid=root.uid).update(reply_count=F("reply_count") - 1)
        messages.success(request, "Moved answer to comment")
        return url

    if action == CLOSE_OFFTOPIC:
        Post.objects.get_all(uid=post.uid).update(status=Post.CLOSED)
        Post.objects.get_all(uid=root.uid).update(reply_count=F("reply_count") - 1)
        content = util.render(name="messages/offtopic_posts.html", user=post.author, comment=comment, posts=post)
        # Create a comment to the post
        Post.objects.create(content=content, type=Post.COMMENT, html=content, parent=post, author=user)
        return url

    if action == DUPLICATE:
        Post.objects.get_all(uid=post.uid).update(status=Post.CLOSED)
        Post.objects.get_all(uid__in=dupes).update(status=Post.CLOSED)
        content = util.render(name="messages/duplicate_posts.html", user=post.author, comment=comment, posts=post)
        # Create a comment to the post
        Post.objects.create(content=content, type=Post.COMMENT, html=content, parent=post, author=user)
        return url

    messages.error(request, "Invalid moderation action given")
    return url


def create_post(title, author, content, post_type, tag_val="", parent=None,root=None, project=None,
                sub_to_root=True):
    "Used to create posts across apps"

    post = Post.objects.create(
        title=title, content=content, tag_val=tag_val,
        author=author, type=post_type, parent=parent, root=root,
        project=project, html=parse_html(content))

    root = root or post.root
    mentioned_users = parse_mentioned_users(content=content)
    subs = Subscription.objects.filter(post=root)

    # Trigger notifications for subscribers and mentioned users
    # async or synchronously
    if tasks.HAS_UWSGI:
        tasks.async_create_sub_messages(subs=subs, author=author, root=root, content=content)
        tasks.async_notify_mentions(users=mentioned_users, root=root, author=author, content=content)
    else:
        tasks.create_sub_messages(subs=subs, author=author, root=root, content=content)
        tasks.notify_mentions(users=mentioned_users, root=root, author=author, content=content)

    # Subscribe the author to the root, if not already
    if sub_to_root:
        create_sub(post=root, sub_type=Subscription.LOCAL_MESSAGE, user=author)

    # Triggers another save in here
    post.add_tags(post.tag_val)

    return post








