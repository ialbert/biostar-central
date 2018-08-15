
import datetime
import logging

from django.utils.timezone import utc
from django.db.models import F
from django.conf import settings
from django.contrib.auth import get_user_model


from .models import Post, Vote, Subscription
from . import util, const

User = get_user_model()


logger = logging.getLogger("engine")


def get_votes(user, thread):

    store = {Vote.BOOKMARK: set(), Vote.UP:set()}

    if user.is_authenticated:
        votes = Vote.objects.filter(post__in=thread, author=user).values_list("post__id", "type")

        for post_id, vote_type in votes:
            store.setdefault(vote_type, set()).add(post_id)

    return store


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


def build_obj_tree(request, obj):

    # Populate the object to build a tree that contains all posts in the thread.
    # Answers sorted before comments.
    user = request.user
    thread = Post.objects.get_thread(obj, user)

    # Build comments tree.
    comment_tree = dict()
    for post in thread:
        if post.is_comment:
            comment_tree.setdefault(post.parent_id, []).append(post)

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

    # Add attributes by mutating the object
    for p in thread:
        decorate(p)
    decorate(obj)

    return comment_tree, thread


def query_topic(user, topic, tag_search=False):
    "Maps known topics to their appropriate querying functions and args."

    if user.is_anonymous:
        tags = ""
    else:
        tags = user.profile.my_tags

    # Acts as a lazy evaluator by only calling a function when its topic is picked
    mapper = {

        const.MYPOSTS: dict(func=Post.objects.my_posts, params=dict(target=user, user=user)),
        const.MYTAGS: dict(func=Post.objects.tag_search, params=dict(text=tags)),
        const.BOOKMARKS: dict(func=Post.objects.my_bookmarks, params=dict(user=user)),
        const.FOLLOWING: dict(func=Post.objects.following, params=dict(user=user)),
        const.LATEST: dict(func=Post.objects.top_level, params=dict(user=user)),

        # "apply" is an added function that can do extra work to the resulting queryset.
        const.VOTES: dict(func=Post.objects.my_post_votes, params=dict(user=user),
                          apply=lambda q: q.distinct()),
        const.OPEN: dict(func=Post.objects.top_level, params=dict(user=user),
                               apply=lambda q: q.filter(type=Post.QUESTION, reply_count=0)),
        const.COMMUNITY: dict(func=User.objects.select_related("profile").all, params=dict()),
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









