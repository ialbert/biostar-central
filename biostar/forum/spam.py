import logging
from django.conf import settings
from django.db.models import Q
from whoosh.writing import AsyncWriter
from whoosh.analysis import StemmingAnalyzer
from whoosh.fields import ID, TEXT, KEYWORD, Schema, NUMERIC
from biostar.forum.models import Post
from biostar.forum import search, auth

logger = logging.getLogger("engine")


def add_to_spam(post, writer):
    writer.update_document(title=post.title,
                           content_length=len(post.content),
                           content=post.content,
                           tags=post.tag_val,
                           uid=post.uid)


def spam_schema():
    analyzer = StemmingAnalyzer()
    schema = Schema(title=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    uid=ID(stored=True),
                    content_length=NUMERIC(stored=True, sortable=True),
                    content=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    tags=KEYWORD(stored=True, commas=True))
    return schema


def report():
    return


def test():
    return


def build_spam_index(overwrite=False):

    # Get all un-indexed spam posts.
    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED), indexed=False)

    # Initialize the spam index
    ix = search.init_index(dirname=settings.SPAM_INDEX_DIR,
                           indexname=settings.SPAM_INDEX_NAME,
                           schema=spam_schema())

    # Index spam posts.
    search.index_posts(posts=spam, ix=ix, overwrite=overwrite, add_func=add_to_spam)

    spam.update(indexed=True)
    logger.info("Built spam index.")
    return ix


def search_spam(post):
    """
    Search spam index for posts similar to this one.
    Returns
    """

    # Initialize the spam index
    ix = search.init_index(dirname=settings.SPAM_INDEX_DIR,
                           indexname=settings.SPAM_INDEX_NAME,
                           schema=spam_schema())

    # Write this post to the index to perform most_like_this()
    writer = AsyncWriter(ix)
    add_to_spam(post=post, writer=writer)
    writer.commit()

    # Search for this post in the spam index
    fields = ['uid']
    results = search.preform_whoosh_search(ix=ix, query=post.uid, fields=fields)

    # Preform more_like_this on this posts content
    similar = results[0].more_like_this('content', top=500)

    # Remove this post from the spam index after results are collected.
    writer = AsyncWriter(ix)
    writer.delete_by_term('uid', text=post.uid)
    writer.commit()

    # Get the results into a list and close the searcher object.
    similar = list(map(search.normalize_result, similar))
    results.searcher.close()

    return similar


def score(post, threshold=None):
    """
    """

    if threshold is None:
        threshold = settings.SPAM_THRESHOLD

    # User's with high enough score automatically given green light.
    #if not post.author.profile.low_rep:
    #    return

    # Search for spam similar to this post.
    similar = search_spam(post=post)

    # Gather the score and calculate the mean.
    scores = [s.score for s in similar]

    if scores:
        mean = sum(scores) / len(scores)
    else:
        mean = 0

    # Update the spam score.
    Post.objects.filter(id=post.id).update(spam_score=mean)

    # If the score exceeds threshold it gets quarantined.
    if mean >= threshold:
        Post.objects.filter(id=post.id).update(spam=Post.SUSPECT)
        auth.log_action(log_text=f"Quarantined post={post.uid}; spam score={mean}")









