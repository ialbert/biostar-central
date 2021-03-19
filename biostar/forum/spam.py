import logging
import os
import shutil
import random
import time
from collections import defaultdict
from math import log, exp
from itertools import groupby, islice, count, chain
from django.conf import settings
from django.db.models import Q
from whoosh.writing import AsyncWriter, BufferedWriter
from whoosh import classify
from whoosh.analysis import StemmingAnalyzer
from whoosh.fields import ID, TEXT, KEYWORD, Schema, NUMERIC, BOOLEAN
from biostar.forum.models import Post
from biostar.accounts.models import Log
from biostar.forum import search, auth, util

logger = logging.getLogger("engine")

HAM_LIMIT = 5000

dirname = os.path.dirname(settings.SPAM_INDEX_DIR)
basename = f"train_{os.path.basename(settings.SPAM_INDEX_DIR)}"
TRAIN_DIR = os.path.join(dirname, basename)

STARTER_UID = 'placeholder'


def sizer(lst, size):
    return size if size <= len(lst) else len(lst)


def spam_schema():
    analyzer = StemmingAnalyzer(cachesize=-1)
    schema = Schema(title=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    content=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    content_length=NUMERIC(stored=True, sortable=True),
                    uid=ID(stored=True),
                    is_spam=BOOLEAN(stored=True))
    return schema


def init_spam_index():
    return search.init_index(dirname=settings.SPAM_INDEX_DIR,
                             indexname=settings.SPAM_INDEX_NAME,
                             schema=spam_schema())


def bootstrap_index(dirname=None, indexname=None):
    """
    Create spam index and add one post from the
    """
    if dirname and indexname:
        ix = search.init_index(dirname=dirname, indexname=indexname, schema=spam_schema())
    else:
        ix = init_spam_index()

    writer = BufferedWriter(ix)
    # Write text to index.
    index_writer(writer=writer, title="Placeholder",
                 content_length=0, is_spam=True,
                 content='CONTENT', uid=STARTER_UID)
    writer.commit()
    writer.close()

    return ix


def index_writer(writer, **kwargs):
    writer.add_document(title=kwargs.get("title"),
                        content_length=kwargs.get("content_length"),
                        content=kwargs.get("content"),
                        uid=kwargs.get("uid"),
                        is_spam=kwargs.get("is_spam", False))


def add_post_to_index(post, writer, is_spam=None):
    """
    Insert post to index.
    """

    index_writer(writer=writer, title=post.title,
                 content_length=len(post.content),
                 content=post.content, is_spam=post.is_spam,
                 uid=f"{post.uid}")


def add_spam(post):

    if post.spam == Post.DEFAULT:
        return

    ix = init_spam_index()
    writer = AsyncWriter(ix)
    add_post_to_index(post=post, writer=writer)
    writer.commit()
    logger.info("Added spam to index.")

    return


def build_spam_index(overwrite=False, add_ham=False, limit=500):
    # Get all un-indexed spam posts.

    spam = Post.objects.filter(spam=Post.SPAM)
    spam = spam.order_by("pk")[:limit]
    spam = list(spam.values_list("id", flat=True))
    # Set indexed flag here so it does not get added to main index.
    if add_ham:
        ham = Post.objects.valid_posts()
        ham = list(ham.values_list("id", flat=True))
        ham = random.sample(ham, k=sizer(ham, size=limit))
    else:
        ham = []

    posts = Post.objects.filter(id__in=chain(spam, ham))

    # Initialize the spam index
    ix = bootstrap_index()

    # Batch index spam posts.
    search.index_posts(posts=posts, ix=ix, overwrite=overwrite, add_func=add_post_to_index)

    logger.info("Built spam index.")
    return ix


def search_spam(post, ix):
    """
    Search spam index for posts similar to this one.
    Returns
    """

    # Add post to index to more_like_this perform search.
    writer = AsyncWriter(ix)
    add_post_to_index(post=post, writer=writer, is_spam=post.is_spam)
    writer.commit()

    # Search for this post in the spam index
    fields = ['uid']
    results = search.perform_search(ix=ix, query=post.uid, fields=fields)

    # Preform more_like_this on this posts content
    similar_content = results[0].more_like_this('content', top=5)

    # Remove this post from the spam index after results are collected.
    writer = AsyncWriter(ix)
    writer.delete_by_term('uid', text=post.uid)
    writer.commit()

    # Get the results into a list and close the searcher object.
    similar_content = list(map(search.copy_hits, similar_content))

    # Close the searcher
    results.searcher.close()

    return similar_content


def remove_spam(post):
    """
    Remove spam from index
    """
    ix = init_spam_index()
    search.remove_post(post=post, ix=ix)
    logger.info("Removed post from spam index.")
    return


def compute_score(post, ix=None):

    ix = ix or init_spam_index()
    N = 1
    bias = -0.25

    # Users above a certain score get green light.
    if not post.author.profile.low_rep:
        return 0

    # Search for spam similar to this post.
    similar_content = search_spam(post=post, ix=ix)

    # Gather the scores for each spam that is similar
    scores = [s.score for s in similar_content if s.is_spam]
    # Take the top N maximum score and compute the mean
    scores = sorted(scores, reverse=True)

    scores = [s for s in scores][:N]

    # Return the mean of the scores.
    if scores:
        value = abs(sum(scores) / len(scores[:N]))
    else:
        # Apply a weighted version of the threshold
        # when no similar posts are found.
        value = settings.SPAM_THRESHOLD + bias

    return value


def score(post, threshold=None):
    """
    """
    if not settings.CLASSIFY_SPAM:
        return

    if threshold is None:
        threshold = settings.SPAM_THRESHOLD

    # User's with high enough score automatically given green light.
    # if not post.author.profile.low_rep:
    #    return

    # Search for spam similar to this post.
    spam_score = compute_score(post=post)

    # Update the spam score.
    Post.objects.filter(id=post.id).update(spam_score=spam_score)

    # If the score exceeds threshold it gets labeled spam.
    if spam_score >= threshold:
        Post.objects.filter(id=post.id).update(spam=Post.SPAM)
        msg = f"auto marked spam :{auth.post_link(post)} spam score={spam_score}"
        auth.db_logger(text=msg)

    return spam_score
