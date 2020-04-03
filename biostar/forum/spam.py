import logging
from math import log
from itertools import groupby, islice, count, chain
from django.conf import settings
from django.db.models import Q
from whoosh.writing import AsyncWriter
from whoosh import classify
from whoosh.analysis import StemmingAnalyzer
from whoosh.fields import ID, TEXT, KEYWORD, Schema, NUMERIC, BOOLEAN
from biostar.forum.models import Post
from biostar.forum import search, auth, util

logger = logging.getLogger("engine")

HAM_LIMIT = 5000

def spam_schema():
    analyzer = StemmingAnalyzer()
    schema = Schema(title=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    content=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    content_length=NUMERIC(stored=True, sortable=True),
                    uid=ID(stored=True),
                    is_spam=BOOLEAN(stored=True))
    return schema


def index_writer(writer, **kwargs):
    writer.update_document(title=kwargs.get("title"),
                           content_length=kwargs.get("content_length"),
                           content=kwargs.get("content"),
                           uid=kwargs.get("uid"),
                           is_spam=kwargs.get("is_spam", False))


def add_post_to_index(post, writer, is_spam=False):
    """
    Insert post to index.
    """
    index_writer(writer=writer, title=post.title,
                 content_length=len(post.content),
                 content=post.content,is_spam=is_spam,
                 uid=post.uid)


def add_file_to_index(fname, delim=","):
    """
    Add file filled with spam and separated by delim.
    """

    # Initialize the index and writer
    ix = search.init_index(dirname=settings.SPAM_INDEX_DIR,
                           indexname=settings.SPAM_INDEX_NAME,
                           schema=spam_schema())
    writer = AsyncWriter(ix)

    # Load the file contents into memory.
    content = open(fname, 'r').read()

    # Split by the delimiter, each item is spam to be indexed.
    content = content.split(delim)
    elapsed, progress = util.timer_func()

    for idx, text in enumerate(content):
        # Print progress
        progress(index=idx, step=500, msg=f"{(idx/len(content) * 100):.02f}%")

        # Write text to index.
        index_writer(writer=writer, title="Placeholder",
                     content_length=len(text), is_spam=True,
                     content=text, uid=util.get_uuid(16))

    # Add some valid posts to the index afrter shuffling with order_by("?")
    ham = Post.objects.filter(Q(spam=Post.NOT_SPAM) | Q(spam=Post.DEFAULT)).order_by("?")[:HAM_LIMIT]
    hcount = ham.count()

    for idx, post in enumerate(ham):
        progress(index=idx, step=500, msg=f"HAM: {(idx/hcount * 100):.02f}%")
        add_post_to_index(post=post, is_spam=False, writer=writer)

    writer.commit()
    elapsed("Committed to index.")


def build_spam_index(overwrite=False):
    # Get all un-indexed spam posts.
    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED), indexed=False)

    # Initialize the spam index
    ix = search.init_index(dirname=settings.SPAM_INDEX_DIR,
                           indexname=settings.SPAM_INDEX_NAME,
                           schema=spam_schema())

    # Index spam posts.
    search.index_posts(posts=spam, ix=ix, overwrite=overwrite, add_func=add_post_to_index)

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
    add_post_to_index(post=post, writer=writer)
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


def compute_score(post, weight=False):

    # Search for spam similar to this post.
    similar = search_spam(post=post)

    # Add weight depending on number of post author has already made.
    # And the user score.
    authored = Post.objects.filter(Q(author=post.author) | Q(lastedit_user=post.author)).count()
    # Scales the values.

    print(post.author.profile.score + authored, authored, post.author.profile.score)
    print(log(1 + (authored + post.author.profile.score) / (len([s.score for s in similar if s.is_spam]) + 5)))
    weight = log(1 + (authored + post.author.profile.score) / (len([s.score for s in similar if s.is_spam]) + 5))

    scores = [s.score - weight for s in similar if s.is_spam]

    # Return the mean of the scores.
    if scores:
        mean = sum(scores) / len(scores)
    else:
        mean = 0
    print(mean)
    return mean


def accuracy(tp, tn, fp, fn):

    return (tp + tn) / (tp + tn + fp + fn)


def specificity(tn, fp):

    return tn / (tn + fp)


def sensitivity(tp, fn):

    return tp / (tp + fn)


def miss_rate(fn, tp):

    return fn / (fn + tp)


def false_discovery_rate(fp, tp):

    return fp / (fp + tp)


def test(threshold=None):

    if threshold is None:
        threshold = settings.SPAM_THRESHOLD

    # Get spam posts.
    spam = Post.objects.filter((Q(spam=Post.SPAM) | Q(status=Post.DELETED)))

    # Get the valid posts and shuffle.
    ham = Post.objects.filter(Q(spam=Post.NOT_SPAM) | Q(status=Post.DEFAULT)).order_by("?")

    ham = ham[:86]

    # Identify spam correctly.
    true_pos = 0

    # Identify valid post correctly.
    true_neg = 0

    # Missed to identify a spam.
    false_neg = 0

    # Mis-identified valid post as spam.
    false_pos = 0

    for post in chain(spam, ham):
        print("-----")
        # Search for spam similar to this post.
        post_score = compute_score(post=post)

        # Positive spam

        if post_score >= threshold:
            if post.is_spam or post.is_deleted:
                #print(post_score, post.is_spam)
                #print(post.content)
                print(post.is_spam or post.is_deleted, "TRUE")
                true_pos += 1
            else:
                #print(post_score, post.is_spam)
                #print(post.content)
                print("FALSE ", post.is_spam or post.is_deleted)

                false_pos += 1
        # Negative spam
        else:
            if post.is_spam or post.is_deleted:
                print("MISSED", post.is_spam or post.is_deleted)
                false_neg += 1
            else:
                true_neg += 1
            pass
        print("-----")
    acc = accuracy(tp=true_pos, tn=true_neg, fp=false_pos, fn=false_neg)
    specif = specificity(tn=true_neg, fp=false_pos)
    sens = sensitivity(tp=true_neg, fn=false_neg)
    missed = miss_rate(fn=false_neg, tp=true_pos)
    false_discovery = false_discovery_rate(fp=false_pos, tp=true_pos)

    print(f"Number of ham: {ham.count()}")
    print(f"Predicted ham: {true_neg}")
    print(f"Number of spam: {spam.count()}")
    print(f"Predicted spam: {true_pos}")

    print(f"Accuracy: {acc}")
    print(f"Specificity : {specif}")
    print(f"Sensitivity : {sens}")
    print(f"Missed rate: {missed}")
    print(f"False discovery rate:{false_discovery}")
    return


def score(post, threshold=None):
    """
    """

    if threshold is None:
        threshold = settings.SPAM_THRESHOLD

    # User's with high enough score automatically given green light.
    # if not post.author.profile.low_rep:
    #    return

    # Search for spam similar to this post.
    post_score = compute_score(post)

    # Update the spam score.
    Post.objects.filter(id=post.id).update(spam_score=post_score)

    # If the score exceeds threshold it gets quarantined.
    if post_score >= threshold:
        Post.objects.filter(id=post.id).update(spam=Post.SUSPECT)
        auth.log_action(log_text=f"Quarantined post={post.uid}; spam score={post_score}")
