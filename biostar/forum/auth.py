
from .models import Post, Tag
from django.contrib import messages

LATEST = "latest"
MYPOSTS, MYTAGS, UNANSWERED, FOLLOWING, BOOKMARKS = "myposts mytags open following bookmarks".split()
POST_TYPES = dict(jobs=Post.JOB, tools=Post.TOOL, tutorials=Post.TUTORIAL,
                  forum=Post.FORUM, planet=Post.BLOG, pages=Post.PAGE)


def post_permissions(request, post):
    """
    Sets permission attributes on a post.
    """
    user = request.user
    is_editable = has_ownership = False

    if user.is_authenticated():

        if user == post.author :
            has_ownership = is_editable = True
        elif user.is_moderator or user.is_staff:
            is_editable = True

    post.is_editable = is_editable
    post.has_ownership = has_ownership

    return post




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



def create_post(title, author, content, tag_val, post_type):

    post = Post.objects.create(
        title=title, content=content, tag_val=tag_val,
        author=author, type=post_type,
    )

    # Triggers a new post save.
    post.add_tags(post.tag_val)

    return post









