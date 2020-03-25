
from biostar.forum.models import Post, Vote
from biostar.forum import search

ACTIVITY_THRESHOLD = 2

SIMILAR_THRESHOLD = 2


def build_spam_index():

    spam = Post.objects.filter(spam=Post.SPAM, indexed=False)

    return


def calc_score(post):

    return


def quarantine(user, post):
    """
    User + post pass a series of tests and rules to determine
    if they might be spam or not
    """

    # User's with high enough score automatically given green light.

    # Mark this post as "maybe" being spam
    Post.objects.filter(id=post.id).update(spam=Post.MAYBE_SPAM)













