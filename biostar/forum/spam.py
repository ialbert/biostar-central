
from biostar.forum.models import Post, Vote
from biostar.forum import search
from biostar.forum.markdown import LINK_PATTERNS

ACTIVITY_THRESHOLD = 2

SIMILAR_THRESHOLD = 2


def build_spam_index():

    spam = Post.objects.filter(spam=Post.SPAM, indexed=False)

    return


def similar_spam(post):



    return


def quarantine(user, post):
    """
    User + post pass a series of tests and rules to determine
    if they might be spam or not
    """

    # User's with high enough score automatically given green light.
    if not user.profile.low_rep:
        return

    # Use has more than the required number of posts.
    post_count = Post.objects.filter(author=user).count()
    vote_count = Vote.objects.filter(author=user).count()
    if (post_count + vote_count) >= ACTIVITY_THRESHOLD:
        return

    # Answer and comments get a green light for now.
    if post.is_answer or post.is_comment:
        return

    # Check the spam index for similar posts.

    # If post exceeds threshold, then make the
    similar = similar_spam(post=post)

    if len(similar) < 5:
        return

    # Mark this post as "maybe" being spam
    Post.objects.filter(id=post.id).update(spam=Post.MAYBE_SPAM)













