
import datetime

from django.contrib import messages
from django.utils.timezone import utc
from django.db.models import F

from .models import Post, Tag, Vote, Subscription


def build_tree(thread, tree={}):

    for post in thread:

        if post.type == Post.COMMENT:
            tree.setdefault(post.parent_id, []).append(post)
    return tree


def get_votes(user, thread):

    store = {Vote.BOOKMARK: set(), Vote.UP:set()}

    if user.is_authenticated:
        pids = [p.id for p in thread]
        votes = Vote.objects.filter(post_id__in=pids, author=user).values_list("post_id", "type")

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



def posts_by_topic(request, topic):
    "Returns a post query that matches a topic"
    user = request.user

    latest = "latest"
    myposts, mytags, unanswered, following, bookmarks, votes = "myposts mytags open following bookmarks votes".split()
    post_types = dict(jobs=Post.JOB, tools=Post.TOOL, tutorials=Post.TUTORIAL,
                      forum=Post.FORUM, planet=Post.BLOG, pages=Post.PAGE)

    # One letter tags are always uppercase
    topic = Tag.fixcase(topic)

    if topic == myposts:
        # Get the posts that the user wrote.
        messages.success(request,'Filtering for posts you contribute to')
        return Post.objects.my_posts(target=user, user=user)

    if topic == mytags:
        # Get the posts that the user wrote.
        messages.success(request,
                         'Posts matching the <b><i class="fa fa-tag"></i> My Tags</b> setting in your user profile')
        return Post.objects.tag_search(user.profile.my_tags)

    if topic == unanswered:
        # Get unanswered posts.
        return Post.objects.top_level(user).filter(type=Post.QUESTION, reply_count=0)

    if topic == following:
        # Get that posts that a user follows.
        messages.success(request, 'Threads that will produce notifications.')
        subs = Subscription.objects.exclude(type=Subscription.NO_MESSAGES).filter(user=user)
        return Post.objects.top_level(user).filter(subs__in=subs)

    if topic == bookmarks:
        # Get that posts that a user bookmarked.
        messages.info(request, f"Showing: {topic}")
        return Post.objects.my_bookmarks(user)

    if topic == votes:
        messages.info(request, f"People voting on your posts")
        return Post.objects.my_post_votes(user).distinct()


    if topic in post_types:
        # A post type.
        return Post.objects.top_level(user).filter(type=post_types[topic])

    if topic and topic != latest:
        # Any type of topic.
        if topic:
            messages.info(request,f"Showing: {topic}" )
        return Post.objects.tag_search(topic)

    # Return latest by default.
    return Post.objects.top_level(user)



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



def create_post(title, author, content, post_type, tag_val="", parent=None,root=None):
    "Used to create posts across apps"


    post = Post.objects.create(
        title=title, content=content, tag_val=tag_val,
        author=author, type=post_type, parent=parent, root=root
    )

    post.add_tags(post.tag_val)

    return post









