from .models import Post, Tag, Vote
from django.contrib import messages


LATEST = "latest"
MYPOSTS, MYTAGS, UNANSWERED, FOLLOWING, BOOKMARKS = "myposts mytags open following bookmarks".split()
POST_TYPES = dict(jobs=Post.JOB, tools=Post.TOOL, tutorials=Post.TUTORIAL,
                  forum=Post.FORUM, planet=Post.BLOG, pages=Post.PAGE)



def build_tree(thread, tree={}):

    for post in thread:
        if post.type == Post.COMMENT:
            tree.setdefault(post.parent_id, []).append(post)
    return


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
    tree = build_tree(thread=thread)
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

    # One letter tags are always uppercase
    topic = Tag.fixcase(topic)

    if topic == MYPOSTS:
        # Get the posts that the user wrote.
        return Post.objects.my_posts(target=user, user=user)

    if topic == MYTAGS:
        # Get the posts that the user wrote.
        messages.success(request,
                         'Posts matching the <b><i class="fa fa-tag"></i> My Tags</b> setting in your user profile')
        return Post.objects.tag_search(user.profile.my_tags)

    if topic == UNANSWERED:
        # Get unanswered posts.
        return Post.objects.top_level(user).filter(type=Post.QUESTION, reply_count=0)

    if topic == FOLLOWING:
        # Get that posts that a user follows.
        messages.success(request, 'Threads that will produce notifications.')
        return Post.objects.top_level(user).filter(subs__user=user)

    if topic == BOOKMARKS:
        # Get that posts that a user bookmarked.
        return Post.objects.my_bookmarks(user)

    if topic in POST_TYPES:
        # A post type.
        return Post.objects.top_level(user).filter(type=POST_TYPES[topic])

    if topic and topic != LATEST:
        # Any type of topic.
        if topic:
            messages.info(request,
                          "Showing: <code>%s</code> &bull; <a href='/'>reset</a>" % topic)
        return Post.objects.tag_search(topic)

    # Return latest by default.
    return Post.objects.top_level(user)



def create_post(title, author, content, post_type, tag_val="", parent=None,root=None):
    "Used to create posts across apps"


    post = Post.objects.create(
        title=title, content=content, tag_val=tag_val,
        author=author, type=post_type, parent=parent, root=root
    )

    post.add_tags(post.tag_val)

    return post









