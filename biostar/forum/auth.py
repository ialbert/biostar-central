import hashlib
import logging
import re
import urllib.parse as urlparse
from datetime import timedelta
from difflib import Differ, SequenceMatcher, HtmlDiff, unified_diff
import bs4
from django.contrib import messages
from django.contrib.auth import get_user_model
from django.core.cache import cache
from django.db import transaction
from django.db.models import F, Q
from django.template import loader
from django.utils.safestring import mark_safe
from django.conf import settings

from biostar.accounts.const import MESSAGE_COUNT
from biostar.accounts.models import Message
from biostar.planet.models import BlogPost, Blog
# Needed for historical reasons.
from biostar.accounts.models import Profile
from biostar.utils.helpers import get_ip
from . import util, awards
from .const import *
from .models import Post, Vote, Subscription, Badge, delete_post_cache, Log, SharedLink, Diff

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


def delete_cache(prefix, user):
    """
    Create key from prefix-user.pk and delete from cache.
    """
    key = f"{prefix}-{user.pk}"

    # Check if it exists and delete object from cache.
    if cache.get(key):
        cache.delete(key)
        logger.debug(f'deleted {key} from cache')

    return


import datetime

ICONS = ["monsterid", "robohash", "wavatar", "retro"]


def gravatar_url(email, style='mp', size=80, force=None):
    global ICONS
    hash_num = hashlib.md5(email).hexdigest()

    # April fools gimmick. Swap icons every hour.
    now = datetime.datetime.now()
    if now.month == 4 and now.day == 1:
        index = now.hour % len(ICONS)
        style = ICONS[index]
        force = True

    data = dict(s=str(size), d=style)

    url = "https://secure.gravatar.com/avatar/%s?" % hash_num

    # Force the new icon.
    if force:
        data['f'] = 'y'

    url += urlparse.urlencode(data)

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
    email = email.encode('utf8', errors="ignore")

    if user.is_anonymous or not user.profile.is_valid:
        # Removes images for suspended users
        email = 'suspended@biostars.org'.encode('utf8')
        style = "monsterid"
        return gravatar_url(email=email, style=style, size=size, force=True)

    # The user has wants a non default icon.
    if user.profile.user_icon != Profile.DEFAULT_ICON:
        style = user.profile.user_icon
        return gravatar_url(email=email, style=style, size=size, force=True)

    # Create the most appropriate default style.
    if user.profile.is_moderator:
        style = "robohash"
    elif user.profile.score > 100:
        style = "retro"
    elif user.profile.score > 0:
        style = "identicon"
    else:
        style = "mp"

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


def create_post(author, title, content, request=None, root=None, parent=None, ptype=Post.QUESTION, tag_val="",
                nodups=True):


    # Check if a post with this exact content already exists.
    post = Post.objects.filter(content=content, author=author).order_by('-creation_date').first()

    # How many seconds since the last post should we disallow duplicates.
    frame = 60
    delta = (util.now() - post.creation_date).seconds if post else frame

    if nodups and delta < frame:
        if request:
            messages.warning(request, "Post with this content was created recently.")
        return post

    post = Post.objects.create(title=title, content=content, root=root, parent=parent,
                               type=ptype, tag_val=tag_val, author=author)

    delete_cache(MYPOSTS, author)
    return post


def diff_ratio(text1, text2):

    # Do not match on spaces
    s = SequenceMatcher(lambda char: re.match(r'\w+', char), text1, text2)
    return round(s.ratio(), 5)


def create_diff(text, post, user):
    """
    Compute and return Diff object for diff between text and post.content
    """

    # Skip on post creation
    if not post:
        return

    ratio = diff_ratio(text1=text, text2=post.content)

    # Skip no changes detected
    if ratio == 1:
        return

    # Compute diff between text and post.
    content = post.content.splitlines()
    text = text.splitlines()

    diff = unified_diff(content, text)
    diff = [f"{line}\n" if not line.endswith('\n') else line for line in diff]
    diff = ''.join(diff)

    # See if a diff has been made by this user in the past 10 minutes
    dobj = Diff.objects.filter(post=post, author=post.author).first()

    # 10 minute time frame between
    frame = 60 * 10
    delta = (util.now() - dobj.created).seconds if dobj else frame

    # Create diff object within time frame or the person editing is a mod.
    if delta >= frame or user != post.author:
        # Create diff object for this user.
        dobj = Diff.objects.create(diff=diff, post=post, author=user)
        post.has_diff = True
        # Only log when anyone but the author commits changes.
        if user != post.author:
            db_logger(user=user, action=Log.EDIT, text=f'edited post', target=post.author, post=post)

    Post.objects.filter(pk=post.pk).update(has_diff=post.has_diff)

    return dobj


def merge_profiles(main, alias):
    """
    Merge alias profile into main
    """

    # Transfer posts
    Post.objects.filter(author=alias).update(author=main)
    Post.objects.filter(lastedit_user=alias).update(lastedit_user=main)

    # Transfer messages
    Message.objects.filter(sender=alias).update(sender=main)
    Message.objects.filter(recipient=alias).update(recipient=main)

    # Do not delete older accounts.
    older = (alias.profile.date_joined < main.profile.date_joined)

    if alias.profile.is_moderator or alias.profile.high_rep or older:
        return

    alias.delete()

    return


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

    # Delete following cache
    delete_cache(FOLLOWING, user)


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
    # Randomly go from one badge to the other
    for award in awards.ALL_AWARDS:

        # Valid award targets the user has earned
        targets = award.get_awards(user)

        for target in targets:
            post = target if isinstance(target, Post) else None
            date = post.lastedit_date if post else user.profile.last_login
            badge = Badge.objects.filter(name=award.name).first()

            valid.append((user, badge, date, post))

    return valid


def get_counts(user):
    # The number of new messages since last visit.
    message_count = Message.objects.filter(recipient=user, unread=True)[:1000].count()

    # The number of new votes since last visit.
    vote_count = Vote.objects.filter(post__author=user, date__gte=user.profile.last_login).exclude(
        author=user)[:1000].count()

    # Planet count since last visit
    planet_count = BlogPost.objects.filter(rank__gte=user.profile.last_login)[:100].count()

    # Spam count since last visit.
    spam_count = Post.objects.filter(spam=Post.SPAM, creation_date__gte=user.profile.last_login)[:1000].count()

    # Moderation actions since last visit.
    mod_count = Log.objects.filter(date__gte=user.profile.last_login)[:100].count()

    # Store the counts into the session.
    counts = dict(mod_count=mod_count, spam_count=spam_count, planet_count=planet_count, message_count=message_count,
                  vote_count=vote_count)

    return counts


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

    # Fetch update the post author score.
    if not post.author == user:
        Profile.objects.filter(user=post.author).update(score=F('score') + change)

    # Calculate counts for the current post
    votes = list(Vote.objects.filter(post=post))
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
        # Reset bookmark cache
        delete_cache(BOOKMARKS, user)

    # Handle accepted vote.
    if vote_type == Vote.ACCEPT:
        Post.objects.filter(uid=post.uid).update(accept_count=accept_count)
        Post.objects.filter(uid=post.root.uid).update(accept_count=F('accept_count') + change)

    return msg, vote, change


def move(request, parent, source, ptype=Post.COMMENT, msg="moved"):
    user = request.user
    url = source.get_absolute_url()

    if source.is_toplevel or not parent:
        return url

    # Move this post to comment of parent
    source.parent = parent
    source.type = ptype

    title = f"{source.get_type_display()}: {source.root.title[:80]}"
    Post.objects.filter(uid=source.uid).update(parent=parent, type=ptype, title=title)

    # Log action and let user know
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg}", post=source)
    source.update_parent_counts()
    return url


def move_post(request, post, parent, **kwargs):
    """
    Move one post to another
    """
    ptype = Post.COMMENT
    msg = f"dragged post to comment"
    return move(request=request,
                parent=parent,
                source=post,
                ptype=ptype,
                msg=msg)


def move_to_answer(request, post, **kwargs):
    """
    Move this post to be an answer
    """

    parent = post.root
    ptype = Post.ANSWER
    msg = "dragged post to answer"
    return move(request=request,
                parent=parent,
                source=post,
                ptype=ptype,
                msg=msg)


def validate_move(user, source, target):
    """
    Return True if moving post from one to another is valid.
    """

    if not source or not target:
        return False

    # cond 1: user is a moderator or author
    valid_user = user.profile.is_moderator or source.author == user

    # cond 2: source and target share the same root
    same_root = source.root.uid == target.root.uid

    # cond 3: source and target are different posts
    is_diff = source.uid != target.uid

    # cond 4: target is not a descendant of source.
    children = set()
    try:
        walk_down_thread(parent=source, collect=children)
        not_desc = (target not in children)
    except Exception as exc:
        logger.error(exc)
        not_desc = False

    # cond 5: source is not top level
    not_toplevel = not source.is_toplevel

    # All conditions need to be met for valid move.
    valid = same_root and is_diff and not_desc and not_toplevel and valid_user

    # Conditions needed to be classified as
    if valid:
        return True

    return False


def db_logger(user, action=Log.MODERATE, text='', target=None, ipaddr=None, post=None):
    """
    Creates a database log.
    """
    Log.objects.create(user=user, action=action, text=text, target=target, ipaddr=ipaddr, post=post)
    logger.info(f"user={user.email} {text} ")
