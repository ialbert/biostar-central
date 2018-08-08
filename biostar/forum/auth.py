
import datetime
import logging

from django.utils.timezone import utc
from django.db.models import F
from django.conf import settings
from django.contrib.auth import get_user_model


from .models import Post, Vote, Subscription, Message
from . import util, const

User = get_user_model()


logger = logging.getLogger("engine")


def build_tree(thread, tree={}):

    for post in thread:

        if post.type == Post.COMMENT:
            tree.setdefault(post.parent_id, []).append(post)
    return tree


def fixcase(name):
    return name.upper() if len(name) == 1 else name.lower()


def get_votes(user, thread):

    store = {Vote.BOOKMARK: set(), Vote.UP:set()}

    if user.is_authenticated:
        pids = [p.id for p in thread]
        votes = Vote.objects.filter(post_id__in=pids, author=user).values_list("post_id", "type")

        for post_id, vote_type in votes:
            store.setdefault(vote_type, set()).add(post_id)

    return store


def get_votes_feed():
    "Returns recent votes feed."

    recent_votes = Vote.objects.filter(type=Vote.UP)[:settings.VOTE_FEED_COUNT]
    # Needs to be put in context of posts
    #recent_votes = post_set.filter(votes__in=recent_votes).order_by("?")


    return


def get_user_feeds():
    "Returns the user location and awards feed."
    return


def post_permissions(request, post):
    """
    Sets permission attributes on a post.
    """
    user = request.user
    is_editable = has_ownership = False

    if user.is_authenticated:

        if user == post.author :
            has_ownership = is_editable = True
        elif user.profile.is_moderator or user.is_staff:
            is_editable = True

    post.is_editable = is_editable
    post.has_ownership = has_ownership

    return post


def build_msg_tree(msg, tree=[]):
    "Build a flat tree with the msg being at the base ( index=0)"

    # Add the current message
    tree.append(msg)

    # Check if it has a parent
    # and recursively add that to that tree.
    if msg.parent_msg and msg.parent_msg != msg:
        tree.append(msg.parent_msg)
        build_msg_tree(msg=msg.parent_msg, tree=tree)

    # End of tree, at the root message
    return tree


def build_obj_tree(request, obj):

    # Populate the object to build a tree that contains all posts in the thread.
    # Answers sorted before comments.
    user = request.user
    thread = [post_permissions(request=request, post=p)
              for p in Post.objects.get_thread(obj, user)]

    # Build tree and gather votes.
    tree = build_tree(thread=thread, tree={})
    votes = get_votes(user=user, thread=thread)

    # Shortcuts to each storage.
    bookmarks = votes[Vote.BOOKMARK]
    upvotes = votes[Vote.UP]

    def decorate(post):
        # Can the current user accept answers
        post.has_bookmark = post.id in bookmarks
        post.has_upvote = post.id in upvotes
        post.can_accept = obj.author == user or post.has_accepted

    # Add attributes by mutating the objects
    map(decorate, thread + [obj])
    # Additional attributes used during rendering
    obj.tree = tree
    obj.answers = [p for p in thread if p.type == Post.ANSWER]

    return obj


def list_message_by_topic(request, topic):

    user = request.user
    # One letter tags are always uppercase
    topic = fixcase(topic)

    if topic == const.MESSAGE:
        return Message.objects.inbox_for(user=user)

    if topic == const.UNREAD:
        return Message.objects.filter(recipient=user, unread=True)

    if topic == const.INBOX:
        return Message.objects.inbox_for(user=user)

    if topic == const.OUTBOX:
        return Message.objects.outbox_for(user=user)

    return Message.objects.inbox_for(user=user)


def list_posts_by_topic(request, topic):
    "Returns a post query that matches a topic"
    user = request.user

    post_types = dict(jobs=Post.JOB, tools=Post.TOOL, tutorials=Post.TUTORIAL,
                      forum=Post.FORUM, planet=Post.BLOG, pages=Post.PAGE)

    # One letter tags are always uppercase
    topic = fixcase(topic)

    if topic == const.MYPOSTS:
        # Get the posts that the user wrote.
        return Post.objects.my_posts(target=user, user=user)

    if topic == const.MYTAGS:
        # Get the posts that the user wrote.
        return Post.objects.tag_search(user.profile.my_tags)

    if topic == const.UNANSWERED:
        # Get unanswered posts.
        return Post.objects.top_level(user).filter(type=Post.QUESTION, reply_count=0)

    if topic == const.FOLLOWING:
        # Get that posts that a user follows.
        subs = Subscription.objects.exclude(type=Subscription.NO_MESSAGES).filter(user=user)
        return Post.objects.top_level(user).filter(subs__in=subs)

    if topic == const.BOOKMARKS:
        # Get that posts that a user bookmarked.
        return Post.objects.my_bookmarks(user)

    if topic == const.VOTES:
        return Post.objects.my_post_votes(user).distinct()

    if topic == const.COMMUNITY:
        # Users that make posts or votes are
        # considered part of the community

        post_set = Post.objects.all()
        users = User.objects.filter(post__in=post_set).distinct()
        return users

    if topic in post_types:
        # A post type.
        return Post.objects.top_level(user).filter(type=post_types[topic])

    if topic and topic != const.LATEST:
        # Any type of topic.
        return Post.objects.tag_search(topic)

    # Return latest by default.
    query = Post.objects.top_level(user)

    return query


def create_sub(post, sub_type, user):
    "Creates a subscription of a user to a post"

    root = post.root
    sub = Subscription.objects.filter(post=root, user=user)
    date = datetime.datetime.utcnow().replace(tzinfo=utc)

    if sub_type == Subscription.DEFAULT_MESSAGES:
        email, local = Subscription.EMAIL_MESSAGE, Subscription.LOCAL_MESSAGE
        sub_type = email if post.is_toplevel else local

    if sub.exists():
        pass
    else:
        sub = Subscription.objects.create(post=root, user=user, type=sub_type, date=date)
        # Increase the subscription count of the root.
        Post.objects.filter(pk=root.pk).update(subs_count=F('subs_count') + 1)

    return sub


def update_vote_count(post):

    vcount = Vote.objects.filter(post=post, type=Vote.UP).count()
    bookcount = Vote.objects.filter(post=post, type=Vote.BOOKMARK).count()

    thread = Post.objects.exclude(status=Post.DELETED).filter(root=post.root, votes__type=Vote.UP)

    thread_score = thread.count()

    # Update the thread score as well
    Post.objects.filter(pk=post.pk).update(vote_count=vcount, book_count=bookcount,
                                           thread_score=thread_score)


def create_messages(body, sender, recipient_list, subject="", parent=None, mtype=None):
    "Create batch message from sender for a given recipient_list"

    subject = subject or f"Message from : {sender.profile.name}"

    msg_list = []
    for rec in recipient_list:

        msg = Message.objects.create(sender=sender, recipient=rec, subject=subject,
                                     body=body, parent_msg=parent, type=mtype)
        msg_list.append(msg)

    return msg_list


def create_vote(author, post, vote_type, updated_type=Vote.EMPTY, update=False):

    vote = Vote.objects.filter(author=author, post=post, type=vote_type).first()

    if update and vote:
        # Update an existing vote type
        Vote.objects.filter(pk=vote.pk).update(type=updated_type)
    else:
        vote = Vote.objects.create(post=post, author=author, type=vote_type)

    # Update the post's vote/bookmark counts
    update_vote_count(post=post)

    return vote


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

    reply_count = json_dict.get("reply_count", 0)
    #thread_score = json_dict.get("thread_score", 0)
    #vote_count = json_dict.get("vote_count", 0)
    view_count = json_dict.get("view_count", 0)

    uid = json_dict.get("id")
    post = Post.objects.filter(uid=uid)
    if post:
        logger.error(f"Post with uid={uid} already exists")
        return post.first()

    post = Post.objects.create(uid=uid, author=author, lastedit_user=lastedit_user,
                               root=root, parent=parent, creation_date=creation_date,
                               lastedit_date=lastedit_date, title=title, has_accepted=has_accepted,
                               type=type, status=status, content=content, html=html, tag_val=tag_val,
                               reply_count=reply_count, #thread_score=thread_score, vote_count=vote_count,
                               view_count=view_count)
    # Trigger another save
    post.add_tags(post.tag_val)

    logger.info(f"Created post.uid={post.uid}")

    return post


def create_post(title, author, content, post_type, tag_val="", parent=None,root=None, project=None):
    "Used to create posts across apps"

    post = Post.objects.create(
        title=title, content=content, tag_val=tag_val,
        author=author, type=post_type, parent=parent, root=root,
        project=project
    )

    # Triggers another save in here
    post.add_tags(post.tag_val)

    return post









