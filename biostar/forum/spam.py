
from biostar.forum.models import Post, Vote
from biostar.forum import search
from biostar.forum.markdown import LINK_PATTERNS

ACTIVITY_THRESHOLD = 2


def has_hyperlink(text):


    return ""


def passes_test(user, post):
    """
    User + post pass a series of tests and rules to determine
    if they might be spam or not
    """

    # User's with high enough score automatically given green light.
    if not user.profile.low_rep:
        return True

    # Use has more than the required number of posts.
    post_count = Post.objects.filter(author=user).count()
    vote_count = Vote.objects.filter(author=user).count()
    if (post_count + vote_count) >= ACTIVITY_THRESHOLD:
        return True

    # User has visited any other page at all.
    
    # Answer and comments get a green light for now.
    if post.is_answer or post.is_comment:
        return True

    if has_hyperlink(post.title):
        return False

    if has_hyperlink(post.content):
        return False

    return False


def is_safe(user, post, request):

    # Test to see if this post is not spam.
    safe = passes_test(user=user, post=post)

    return safe





