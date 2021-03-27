import hashlib
import logging
import urllib.parse as urlparse

from django.contrib import messages
from django.contrib.auth import get_user_model
from django.core.cache import cache
from django.db import transaction
from django.db.models import F, Q
from django.template import loader
from django.utils.safestring import mark_safe

from biostar.accounts.const import MESSAGE_COUNT
from biostar.accounts.models import Message
# Needed for historical reasons.
from biostar.accounts.models import Profile
from . import util, awards
from .const import *
from .models import Post, Vote, Subscription, Badge, delete_post_cache, Log

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


def gravatar_url(email, style='mp', size=80, force=None):
    hash_num = hashlib.md5(email).hexdigest()

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


def create_post(author, title, content, root=None, parent=None, ptype=Post.QUESTION, tag_val=""):
    # Check if a post with this exact content already exists.
    post = Post.objects.filter(content=content, author=author, is_toplevel=True).first()
    if post:
        logger.info("Post with this content already exists.")
        return post

    post = Post.objects.create(title=title, content=content, root=root, parent=parent,
                               type=ptype, tag_val=tag_val, author=author)

    delete_cache(MYPOSTS, author)
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
    planet_count = 0

    # Spam count since last visit.
    spam_count = 0

    # Moderation actions since last visit.
    mod_count = 0

    # Store the counts into the session.
    counts = dict(mod_count=mod_count, spam_count=spam_count, planet_count=planet_count)

    # TODO: needs to be changed
    counts[MESSAGE_COUNT] = message_count
    counts[VOTES_COUNT] = vote_count

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
        # Reset bookmark cache
        delete_cache(BOOKMARKS, user)

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


def removal_condition(post, user, age=1):
    """
    Removal condition for the post.
    """

    # Only authors may remove their own posts
    if post.author != user:
        return False

    # Post older than a day may not be removed
    if post.age_in_days > age:
        return False

    # If the post has children it may not be removed
    if Post.objects.filter(parent=post).exclude(pk=post.id):
        return False

    # If the post has votes it may not be removed
    if post.vote_count:
        return False

    return True


def delete_post(request, post, **kwargs):
    """
    Post may be marked as deleted or removed entirely
    """
    user = request.user

    # Decide on the removal
    remove = removal_condition(post, user)

    if remove:
        msg = f"removed post"
        messages.info(request, mark_safe(msg))
        db_logger(user=user, post=post, text=msg)
        post.delete()
        url = "/"
    else:
        Post.objects.filter(uid=post.uid).update(status=Post.DELETED)
        post.recompute_scores()
        msg = f"deleted post"
        messages.info(request, mark_safe(msg))
        db_logger(user=user, post=post, text=msg)
        url = post.get_absolute_url()

    # Recompute post score.
    if not post.is_toplevel:
        post.root.recompute_scores()

    return url


def open(request, post, **kwargs):
    if post.is_spam and post.author.profile.low_rep:
        post.author.profile.bump_over_threshold()

    user = request.user
    Post.objects.filter(uid=post.uid).update(status=Post.OPEN, spam=Post.NOT_SPAM)
    post.recompute_scores()

    post.root.recompute_scores()
    msg = f"opened post"
    url = post.get_absolute_url()
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg}", post=post)
    return url


def bump(request, post, **kwargs):
    now = util.now()
    user = request.user

    Post.objects.filter(uid=post.uid).update(lastedit_date=now, rank=now.timestamp())
    msg = f"bumped post"
    url = post.get_absolute_url()
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg}", post=post)

    return url


def change_user_state(mod, target, state):
    """
    Changes user state.
    """

    # Only moderators may change user states.
    if not mod.profile.is_moderator:
        logger.error(f"{mod} is not a moderator")
        return

    # Cannot moderate self.
    if mod == target:
        logger.error(f"{mod} cannot moderate themselves")
        return

    # The target may not be a moderator.
    if target.profile.is_moderator:
        logger.info(f"{mod} cannot alter state on a moderator {target}")
        return

    # Set the new state.
    target.profile.state = state
    target.profile.save()

    # Generate the logging message.
    msg = f"changed user state to {target.profile.get_state_display()}"
    db_logger(user=mod, action=Log.MODERATE, target=target, text=msg, post=None)


def toggle_spam(request, post, **kwargs):
    """
    Toggles spam status on post based on a status
    """

    url = post.get_absolute_url()

    # Moderators may only be suspended by admins (TODO).
    if post.author.profile.is_moderator and post.spam in (Post.DEFAULT, Post.NOT_SPAM):
        messages.warning(request, "cannot toggle spam on a post created by a moderator")
        return url

    # The user performing the action.
    user = request.user

    # Drop the cache for the post.
    delete_post_cache(post)

    # Current state of the toggle.
    if post.is_spam:
        Post.objects.filter(id=post.id).update(spam=Post.NOT_SPAM, status=Post.OPEN)
    else:
        Post.objects.filter(id=post.id).update(spam=Post.SPAM, status=Post.CLOSED)

    # Refetch up to date state of the post.
    post = Post.objects.filter(id=post.id).get()

    # Set the state for the user (only non moderators are affected)
    state = Profile.SUSPENDED if post.is_spam else Profile.NEW

    # Apply the user change.
    change_user_state(mod=user, target=post.author, state=state)

    # Generate logging messages.
    if post.is_spam:
        text = f"marked post as spam"
    else:
        text = f"restored post from spam"

        # Set indexed flag to False, so it's removed from spam index
        Post.objects.filter(id=post.id).update(indexed=False)

    # Set a logging message.
    messages.success(request, text)

    # Submit the log into the database.
    db_logger(user=user, action=Log.MODERATE, target=post.author, text=text, post=post)

    url = post.get_absolute_url()

    return url


def move(request, parent, source, ptype=Post.COMMENT, msg="moved"):
    user = request.user
    url = source.get_absolute_url()

    if source.is_toplevel or not parent:
        return url

    # Move this post to comment of parent
    source.parent = parent
    source.type = ptype
    source.save()

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
    msg = f"moved post"
    return move(request=request,
                parent=parent,
                source=post,
                ptype=ptype,
                msg=msg)


def move_to_comment(request, post, **kwargs):
    """
    move this post to a comment
    """

    parent = post.root
    ptype = Post.COMMENT
    msg = "moved comment"
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
    msg = "moved answer"
    return move(request=request,
                parent=parent,
                source=post,
                ptype=ptype,
                msg=msg)


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
    msg = "closed"
    url = rationale.get_absolute_url()
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg}", post=post)
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
    msg = "duplicate"
    messages.info(request, mark_safe(msg))
    db_logger(user=user, text=f"{msg}", post=post)
    return url


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


def off_topic(request, post, **kwargs):
    """
    Marks post as off topic. Generate off topic comment.
    """
    user = request.user
    if "offtopic" not in post.tag_val:
        post.tag_val += ",offtopic "
        post.save()

        # Generate off comment.
        content = "This post is off topic."
        create_post(ptype=Post.COMMENT, parent=post, content=content, title='', author=request.user)
        msg = "off topic"
        messages.info(request, mark_safe(msg))
        db_logger(user=user, text=f"{msg}", post=post)
    else:
        messages.warning(request, "post has been already tagged as off topic")

    url = post.get_absolute_url()
    return url


def relocate(request, post, **kwds):
    """
    Moves an answer to a comment and viceversa.
    """
    url = post.get_absolute_url()

    if post.is_toplevel:
        messages.warning(request, "cannot relocate a top level post")
        return url

    if post.type == Post.COMMENT:
        msg = f"moved comment to answer"
        post.type = Post.ANSWER
    else:
        msg = f"moved answer to comment"
        post.type = Post.COMMENT

    post.parent = post.root
    post.save()

    db_logger(user=request.user, post=post, text=f"{msg}")
    messages.info(request, msg)
    return url


def moderate(request, post, action, parent, comment=""):
    # Bind an action to a function.
    action_map = {
        REPORT_SPAM: toggle_spam,
        DUPLICATE: duplicated,
        BUMP_POST: bump,
        OPEN_POST: open,
        OFF_TOPIC: off_topic,
        DELETE: delete_post,
        CLOSE: close,
        MOVE_COMMENT: move_to_comment,
        MOVE_ANSWER: move_to_answer,

        # TODO: change to words!
        100: relocate,
    }

    if action in action_map:
        mod_func = action_map[action]
        url = mod_func(request=request, post=post, parent=parent, comment=comment)
    else:
        url = post.get_absolute_url()
        msg = "Unknown moderation action given."
        logger.error(msg)

    return url


def db_logger(user=None, action=Log.MODERATE, text='', target=None, ipaddr=None, post=None):
    """
    Creates a database log.
    """
    Log.objects.create(user=user, action=action, text=text, target=target, ipaddr=ipaddr, post=post)
    logger.info(f"user={user.email} {text} ")
