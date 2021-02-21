import datetime
import logging
import json
import hashlib

import urllib.parse as urlparse
from urllib import request

from django.template import loader
from django.conf import settings
from django.contrib import messages
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import F, Q
from django.utils.timezone import utc
from django.core.cache import cache
from django.core.paginator import Paginator
from django.shortcuts import reverse
from biostar.accounts.models import Profile
from . import util, awards
from .const import *
from .models import Post, Vote, PostView, Subscription, Award, Badge
from django.utils.safestring import mark_safe


# Needed for historical reasons.
from biostar.accounts.auth import db_logger

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


def convert_html():
    """
    Converts html to text
    """
    return


def gravatar_url(email, style='mp', size=80):
    hash_num = hashlib.md5(email).hexdigest()

    url = "https://secure.gravatar.com/avatar/%s?" % hash_num
    url += urlparse.urlencode({
        's': str(size),
        'd': style,
    }
    )
    return url


def encode_email(email, key):
    """
    Use key to encode email
    """
    return


def decode_email(email):
    """
    Use api key to decode email
    """
    return


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


def create_post_from_json(**json_data):

    post_uid = json_data['id']

    # Check to see if the uid already exists
    post = Post.objects.filter(uid=post_uid).first()

    # Update an existing post
    if post:

        post.content = json_data['']
        post.lastedit_date = json_data['lastedit_date']
        post.creation_date = json_data['creation_date']

    # data = {
    #     'id': self.id,
    #     'uid': self.uid,
    #     'title': self.title,
    #     'type': self.get_type_display(),
    #     'type_id': self.type,
    #     'creation_date': util.datetime_to_iso(self.creation_date),
    #     'lastedit_date': util.datetime_to_iso(self.lastedit_date),
    #     'lastedit_user_id': self.lastedit_user.id,
    #     'author_id': self.author.id,
    #     'author_uid': self.author.profile.uid,
    #     'lastedit_user_uid': self.lastedit_user.profile.uid,
    #     'author': self.author.name,
    #     'status': self.get_status_display(),
    #     'status_id': self.status,
    #     'thread_score': self.thread_votecount,
    #     'rank': self.rank,
    #     'vote_count': self.vote_count,
    #     'view_count': self.view_count,
    #     'reply_count': self.reply_count,
    #     'comment_count': self.comment_count,
    #     'book_count': self.book_count,
    #     'subs_count': self.subs_count,
    #     'answer_count': self.root.reply_count,
    #     'has_accepted': self.has_accepted,
    #     'parent_id': self.parent.id,
    #     'root_id': self.root_id,
    #     'xhtml': self.html,
    #     'content': self.content,
    #     'tag_val': self.tag_val,
    #     'url': f'{settings.PROTOCOL}://{settings.SITE_DOMAIN}{self.get_absolute_url()}',
    # }
    
    return


def create_post(author, title, content, root=None, parent=None, ptype=Post.QUESTION, tag_val=""):

    # Check if a post with this content already exists.
    post = Post.objects.filter(content=content, author=author).first()
    if post:
        logger.info("Post with this content already exists.")
        return post

    post = Post.objects.create(title=title, content=content, root=root, parent=parent,
                               type=ptype, tag_val=tag_val, author=author)

    return post


def create_subscription(post, user, sub_type=None, update=False):
    """
    Creates subscription to a post. Returns a list of subscriptions.
    """
    subs = Subscription.objects.filter(post=post.root, user=user)
    sub = subs.first()

    default = Subscription.TYPE_MAP.get(user.profile.message_prefs,
                                        Subscription.LOCAL_MESSAGE)

    empty = sub_type is None
    # Get the current sub type from what's given or the existing sub
    sub_type = None if empty else sub_type
    # No type has been given so default
    sub_type = sub_type or default

    # Ensure the sub type is not set to something wrote
    if sub and update:
        # Update an existing subscription
        sub.type = sub_type
        sub.save()
    else:
        # Drop all existing subscriptions for the user by default.
        subs.delete()
        Subscription.objects.create(post=post.root, user=user, type=sub_type)

    # Recompute subscription count
    subs_count = Subscription.objects.filter(post=post.root).exclude(type=Subscription.NO_MESSAGES).count()

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
    query = Post.objects.valid_posts(u=user, root=root).exclude(pk=root.id)

    # Filter spam/deleted comments or answers.
    if user.is_anonymous or not user.profile.is_moderator:
        query = query.exclude(Q(status=Post.DELETED) | Q(spam=Post.SPAM))

    query = query.select_related("lastedit_user__profile", "author__profile", "root__author__profile")

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
        if user.is_authenticated:
            post.can_accept = not post.is_toplevel and (user == post.root.author or user.profile.is_moderator)
            post.can_moderate = user.profile.is_moderator
            post.is_editable = (user == post.author or user.profile.is_moderator)
        else:
            post.can_accept = False
            post.is_editable = False
            post.can_moderate = False

        return post

    # Decorate the objects for easier access
    thread = list(map(decorate, thread))

    # Decorate the root post
    root = decorate(root)

    # Select the answers from the thread.
    answers = [p for p in thread if p.type == Post.ANSWER]

    return root, comment_tree, answers, thread


def valid_awards(user):
    """
    Return list of valid awards for a given user
    """

    valid = []
    for award in awards.ALL_AWARDS:

        # Valid award targets the user has earned
        targets = award.validate(user)
        for target in targets:

            date = util.now()
            post = target if isinstance(target, Post) else None
            badge = Badge.objects.filter(name=award.name).first()

            # Do not award a post multiple times.
            already_awarded = Award.objects.filter(user=user, badge=badge, post=post).exists()
            if post and already_awarded:
                continue

            valid.append((user, badge, date, post))

    return valid


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
        msg = f"{vote.get_type_display()} removed"
        change = -1
        vote.delete()
    else:
        change = +1
        vote = Vote.objects.create(author=user, post=post, type=vote_type)
        msg = f"{vote.get_type_display()} added"

    if post.author == user:
        # Author making the change
        change = 0
        return msg, vote, change

    # Fetch update the user score.
    Profile.objects.filter(user=post.author).update(score=F('score') + change)

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


def mod_rationale(post, user, template, ptype=Post.ANSWER, extra_context=dict()):
    tmpl = loader.get_template(template)
    context = dict(user=post.author)
    context.update(extra_context)
    content = tmpl.render(context)

    # Load answer explaining post being off topic.
    post = Post.objects.create(content=content, type=ptype, parent=post, root=post.root, author=user)

    return post


def only_delete(post, user):
    """
    Return if this post should get a 'Delete' flag instead of being removed from database.
    """
    children = Post.objects.filter(parent_id=post.id).exclude(pk=post.id)
    # The conditions where post can only be deleted.
    cond1 = children or post.age_in_days > 7
    cond2 = post.vote_count > 1 or (post.author != user)

    only_del = cond1 or cond2

    return only_del


def delete_post(post, user):

    if only_delete(post, user):
        # Deleted posts can be un=deleted by re-opening them.
        Post.objects.filter(uid=post.uid).update(status=Post.DELETED)
        url = post.root.get_absolute_url()
        msg = f"Deleted post: {post_link(post)}"
        post.recompute_scores()
        post.root.recompute_scores()
        return url, msg

    # Redirect depends on the level of the post.
    if post.is_toplevel:
        url = "/"
    else:
        url = post.root.get_absolute_url()
        post.root.recompute_scores()

    # Remove post from the database with no trace.
    msg = f"Removed post: {post_link(post)}"
    post.delete()

    return url, msg


def open(request, post, **kwargs):

    if post.is_spam and post.author.profile.low_rep:
        post.author.profile.bump_over_threshold()

    user = request.user
    Post.objects.filter(uid=post.uid).update(status=Post.OPEN, spam=Post.NOT_SPAM)
    post.recompute_scores()

    post.root.recompute_scores()
    msg = f"Opened post: {post_link(post)}"
    url = post.get_absolute_url()
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg} ; post.uid={post.uid}.")
    return url


def bump(request, post, **kwargs):
    now = util.now()
    user = request.user

    Post.objects.filter(uid=post.uid).update(lastedit_date=now, rank=now.timestamp())
    msg = f"Post bumped : {post_link(post)}"
    url = post.get_absolute_url()
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg} ; post.uid={post.uid}.")

    return url


def post_link(post):
    url = post.get_absolute_url()
    link = f'<a href="{url}">{post.title} ({post.uid})</a>'
    return link


def user_link(user):
    url = user.get_absolute_url()
    link = f'<a href="{url}">{user.name} ({user.uid})</a>'
    return link


def toggle_spam(request, post, **kwargs):
    """
    Toggles spam status on post based on a status
    """

    # The user performing the action.
    user = request.user

    # The spam status set by the toggle.
    spam = Post.NOT_SPAM if post.is_spam else Post.SPAM

    # Update the object bypassing the signals.
    post.spam = spam
    Post.objects.filter(id=post.id).update(spam=spam)

    # Moderators may only be suspended by admins (TODO).
    if post.author.profile.is_moderator:
        messages.warning(request, "cannot toggle spam on a post created by a moderator")
        return

    # Set the state for the user.
    post.author.profile.state = Profile.SUSPENDED if post.is_spam else Profile.NEW
    post.author.profile.save()

    # Generate logging messages.
    if post.is_spam:
        text = f'marked post {post_link(post)} as spam'
    else:
        text = f'restored post {post_link(post)} from spam'

    # Submit the log into the database.
    db_logger(user=user, text=text)
    messages.info(request, mark_safe(text))
    url = post.get_absolute_url()
    return url


def close(request, post, comment, **kwargs):
    """
    Close this post and provide a rationale for closing as well.
    """
    user = request.user
    Post.objects.filter(uid=post.uid).update(status=Post.CLOSED)
    # Generate a rationale post on why this post is closed.
    context = dict(comment=comment)
    rationale = mod_rationale(post=post, user=user,
                              template="messages/closed.md",
                              extra_context=context)
    msg = f"Closed {post_link(post)}. "
    url = rationale.get_absolute_url()
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg} ; post.uid={post.uid}.")
    return url


def duplicated(request, post, comment, **kwargs):

    # Generate a rationale post on why this post is a duplicate.
    Post.objects.filter(uid=post.uid).update(status=Post.CLOSED)
    user = request.user

    dupes = comment.split("\n")[:5]
    dupes = list(filter(lambda d: len(d), dupes))
    context = dict(dupes=dupes, comment=comment)
    rationale = mod_rationale(post=post, user=user,
                              template="messages/duplicate_posts.md",
                              extra_context=context)
    url = rationale.get_absolute_url()
    msg = f"Closed duplicated post {post_link(post)}."
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg} ; post.uid={post.uid}.")
    return url


def delete(request, post, **kwargs):
    """
    Delete this post or complete remove it from the database.
    """
    user = request.user
    url, msg = delete_post(post=post, user=user)
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg} ; post.uid={post.uid}.")

    return url


def moderate(request, post, action, comment=""):

    # Bind an action to a function.
    action_map = {REPORT_SPAM: toggle_spam,
                  DUPLICATE: duplicated,
                  BUMP_POST: bump,
                  OPEN_POST: open,
                  DELETE: delete,
                  CLOSE: close}

    if action in action_map:
        mod_func = action_map[action]
        url = mod_func(request=request, post=post, comment=comment)
    else:
        url = post.get_absolute_url()
        msg = "Unknown moderation action given."
        logger.error(msg)

    return url
