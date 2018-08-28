
import bleach
import datetime
import logging
import re
import mistune

from django.contrib import messages
from django.utils.timezone import utc
from django.db.models import F
from django.conf import settings
from django.contrib.auth import get_user_model
from biostar.message import tasks
from biostar.utils.shortcuts import reverse
from .models import Post, Vote, Subscription
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



def build_obj_tree(request, obj):

    # Populate the object to build a tree that contains all posts in the thread.
    # Answers sorted before comments.
    user = request.user
    thread = Post.objects.get_thread(obj, user)

    # Gather votes
    votes = get_votes(user=user, thread=thread)

    # Shortcuts to each storage.
    bookmarks = votes[Vote.BOOKMARK]
    upvotes = votes[Vote.UP]

    def decorate(post):
        # Can the current user accept answers
        post.has_bookmark = post.id in bookmarks
        post.has_upvote = post.id in upvotes
        post.can_accept = obj.author == user or post.has_accepted

    # Build comments tree.
    comment_tree = dict()
    decorate(obj)

    for post in thread:
        if post.is_comment:
            comment_tree.setdefault(post.parent_id, []).append(post)
        # Add attributes by mutating the object
        decorate(post=post)

    return comment_tree, thread


def query_topic(user, topic, tag_search=False):
    "Maps known topics to their appropriate querying functions and args."

    if user.is_anonymous:
        tags = ""
    else:
        tags = user.profile.my_tags

    # Acts as a lazy evaluator by only calling a function when its topic is picked
    mapper = {

        MYPOSTS: dict(func=Post.objects.my_posts, params=dict(target=user, user=user)),
        MYTAGS: dict(func=Post.objects.tag_search, params=dict(text=tags)),
        BOOKMARKS: dict(func=Post.objects.my_bookmarks, params=dict(user=user)),
        FOLLOWING: dict(func=Post.objects.following, params=dict(user=user)),
        LATEST: dict(func=Post.objects.top_level, params=dict(user=user)),

        # "apply" is an added function that can do extra work to the resulting queryset.
        VOTES: dict(func=Post.objects.my_post_votes, params=dict(user=user),
                          apply=lambda q: q.distinct()),
        OPEN: dict(func=Post.objects.top_level, params=dict(user=user),
                               apply=lambda q: q.filter(type=Post.QUESTION, reply_count=0)),
        COMMUNITY: dict(func=User.objects.select_related("profile").all, params=dict()),
    }

    if mapper.get(topic):
        func, params = mapper[topic]["func"], mapper[topic].get("params")
        apply_extra = mapper[topic].get("apply", lambda q: q)
        query = apply_extra(func(**params))
    else:
        query = None

    # Query any topic as a tag
    if tag_search and query is None:
        query = Post.objects.tag_search(topic)

    return query


def list_posts_by_topic(request, topic):
    "Returns a post query that matches a topic"
    user = request.user

    post_types = dict(jobs=Post.JOB, tools=Post.TOOL, tutorials=Post.TUTORIAL,
                      forum=Post.FORUM, planet=Post.BLOG, pages=Post.PAGE)

    # One letter tags are always uppercase.
    topic = util.fixcase(topic)
    query = query_topic(user=user, topic=topic, tag_search=True)

    # A post type.
    if topic in post_types:
        query = Post.objects.top_level(user).filter(type=post_types[topic])

    # Return latest by default.
    if query is None:
        query = Post.objects.top_level(user)

    return query


def create_sub(post,  user, sub_type=None):
    "Creates a subscription of a user to a post"

    root = post.root
    sub = Subscription.objects.filter(post=root, user=user)
    date = datetime.datetime.utcnow().replace(tzinfo=utc)
    exists = sub.exists()

    # Subscription already exists
    if exists and sub_type is None:
        return sub
    # Update an existing object
    elif exists:
        Subscription.objects.update(type=sub_type)
        # The sub is being changed to "No message"
        if sub_type == Subscription.NO_MESSAGES:
            Post.objects.filter(pk=root.pk).update(subs_count=F('subs_count') - 1)

    # Create a new object
    else:
        sub = Subscription.objects.create(post=root, user=user, type=sub_type, date=date)
        # Increase the subscription count of the root.
        Post.objects.filter(pk=root.pk).update(subs_count=F('subs_count') + 1)

    return sub


def update_vote_count(post):

    vcount = Vote.objects.filter(post=post, type=Vote.UP).count()
    bookcount = Vote.objects.filter(post=post, type=Vote.BOOKMARK).count()
    # TODO:Thread score compted wrong here
    thread = Post.objects.exclude(status=Post.DELETED).filter(root=post.root, votes__type=Vote.UP)

    thread_score = thread.count()

    # Update the thread score as well
    Post.objects.filter(pk=post.pk).update(vote_count=vcount, book_count=bookcount,
                                           thread_score=thread_score)


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
    if post.exists() or status == Post.DELETED:
        logger.error(f"Post with uid={uid} already exists or status is deleted.")
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
        url = reverse("public_profile", kwargs=dict(uid=user.profile.uid))
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
    post.update_reply_count()

    return url


def moderate_post(request, action, post, comment=None, dupes=[]):
    """Used to moderate a post given a specific action"""

    root = post.root
    user = request.user
    now = datetime.datetime.utcnow().replace(tzinfo=utc)
    url = post.root.get_absolute_url()
    query_func = Post.objects.get_all

    root_has_accepted = lambda: query_func(root=root, type=Post.ANSWER, has_accepted=True).count()

    # Acts as a lazy evaluator by only calling a function when its action is picked
    action_map = {

        BUMP_POST: dict(func=query_func(uid=post.uid).update,
                        params=dict(lastedit_date=now, lastedit_user=request.user),
                        msg="Post bumped"),

        MOD_OPEN: dict(func=query_func(uid=post.uid).update,
                   params=dict(status=Post.OPEN),
                   msg=f"Opened post: {post.title}"),

        DELETE: dict(func=delete_post, params=dict(post=post, request=request)),

        CROSSPOST: dict(content="messages/crossposted.html"),

        TOGGLE_ACCEPT: dict(func=query_func(uid=post.uid).update,
                            params=dict(has_accepted=not post.has_accepted),
                            apply=lambda q: query_func(uid=root.uid).update(has_accepted=root_has_accepted())),

        MOVE_TO_ANSWER: dict(func=query_func(uid=post.uid).update,
                             params=dict(type=Post.ANSWER, parent=post.root),
                             msg="Moved comment to answer",
                             apply=lambda q: query_func(uid=root.uid).update(reply_count=F("reply_count") + 1)),

        MOVE_TO_COMMENT: dict(func=query_func(uid=post.uid).update,
                              params=dict(type=Post.COMMENT, parent=post.root),
                              msg="Moved answer to comment",
                              apply=lambda q: query_func(uid=root.uid).update(reply_count=F("reply_count") - 1)),

        CLOSE_OFFTOPIC: dict(func=query_func(uid=post.uid).update,
                             params=dict(status=Post.CLOSED),
                             msg=f"Closed post: {post.title}",
                             content = "messages/offtopic_posts.html",
                             apply=lambda q: query_func(uid=root.uid).update(reply_count=F("reply_count") - 1)),

        DUPLICATE: dict(func=query_func(uid=post.uid).update,
                        params=dict(status=Post.CLOSED),
                        apply=lambda q: query_func(uid=post.uid).filter(uid__in=dupes),
                        content="messages/duplicate_posts.html"),
    }

    # Valid moderation action is being taken

    if action_map.get(action):
        # Specific action to moderate
        action = action_map[action].get

        # Function and params associated with it
        func = action("func", lambda x: x)
        params = action("params", dict(x=None))
        msg, extra_content = action("msg"), action("content")
        extra_protocol = action("apply", lambda x: x)
        # Apply the moderation action
        output = extra_protocol(func(**params))

        # Get the correct url to redirect to ( matters most when a root post is deleted).
        url = output if action == DELETE else url

        if msg:
            messages.success(request, msg)
        if extra_content:
            content = util.render(name=extra_content, user=post.author, comment=comment, posts=output or post)
            # Create a comment to the post
            Post.objects.create(content=content, type=Post.COMMENT, html=content, parent=post, author=user)
    else:
        messages.error(request, "Invalid moderation action given")

    return url


def create_post(title, author, content, post_type, tag_val="", parent=None,root=None, project=None,
                sub_to_root=True):
    "Used to create posts across apps"

    post = Post.objects.create(
        title=title, content=content, tag_val=tag_val,
        author=author, type=post_type, parent=parent, root=root,
        project=project, html=parse_html(content),
    )

    root = root or post.root
    # Trigger notifications for subscribers and mentioned users
    # async or synchronously

    mentioned_users = parse_mentioned_users(content=content)
    mention_params = dict(users=mentioned_users, root=root, author=author, content=content)
    subs = Subscription.objects.filter(post=root)
    subs_params = dict(subs=subs, author=author, root=root, content=content)

    if tasks.HAS_UWSGI:
        tasks.async_create_sub_messages(**subs_params)
        tasks.async_notify_mentions(**mention_params)
    else:
        tasks.create_sub_messages(**subs_params)
        tasks.notify_mentions(**mention_params)

    # Subscribe the author to the root, if not already
    if sub_to_root:
        create_sub(post=root, sub_type=Subscription.LOCAL_MESSAGE, user=author)

    # Triggers another save in here
    post.add_tags(post.tag_val)

    return post








