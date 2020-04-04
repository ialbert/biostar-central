import logging
import random
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


STARTER_UID = 'placeholder'


def spam_schema():
    analyzer = StemmingAnalyzer()
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


def bootstrap_index():
    """
    Create spam index and add one post from the
    """
    ix = init_spam_index()
    writer = AsyncWriter(ix)
    # Write text to index.
    index_writer(writer=writer, title="Placeholder",
                 content_length=0, is_spam=True,
                 content='', uid=STARTER_UID)

    return


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


def add_file_to_index(fname, delim=",", train_test_split=.2):
    """
    Add file filled with spam and separated by delim.
    """

    # Initialize the index and writer
    ix = init_spam_index()
    writer = AsyncWriter(ix)

    bootstrap_index()

    # Load the file contents into memory.
    content = open(fname, 'r').read()

    # Split by the delimiter, each item is spam to be indexed.
    content = content.split(delim)
    split = int(len(content) * train_test_split)
    random.shuffle(content)
    train = content[split:]
    test = content[:split]

    elapsed, progress = util.timer_func()

    for idx, text in enumerate(train):
        # Print progress
        progress(index=idx, step=500, msg=f"{(idx/len(train) * 100):.02f}%")

        # Write text to index.
        index_writer(writer=writer, title="Placeholder",
                     content_length=len(text), is_spam=True,
                     content=text, uid=util.get_uuid(16))

    # Add some valid posts to the index afrter shuffling with order_by("?")
    ham = Post.objects.valid_posts().order_by("?")[:HAM_LIMIT]
    hcount = ham.count()

    for idx, post in enumerate(ham):
        progress(index=idx, step=500, msg=f"HAM: {(idx/hcount * 100):.02f}%")
        add_post_to_index(post=post, is_spam=False, writer=writer)

    writer.commit()
    elapsed("Committed to index.")
    return test


def build_spam_index(overwrite=False):
    # Get all un-indexed spam posts.
    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED), indexed=False)

    # Initialize the spam index
    ix = init_spam_index()

    bootstrap_index()

    # Index spam posts.
    search.index_posts(posts=spam, ix=ix, overwrite=overwrite, add_func=add_post_to_index)

    spam.update(indexed=True)
    logger.info("Built spam index.")
    return ix


def search_spam(text):
    """
    Search spam index for posts similar to this one.
    Returns
    """

    # Initialize the spam index
    ix = init_spam_index()

    # Preform more like this on the "TEXT"

    # Search for this post in the spam index
    results = search.preform_whoosh_search(ix=ix, query=STARTER_UID, fields=['uid'])

    # Preform more_like_this on this posts content
    similar = results[0].more_like_this(text=text, top=500)

    # Get the results into a list and close the searcher object.
    similar = list(map(search.normalize_result, similar))
    results.searcher.close()

    return similar


def compute_score(text, post=None):

    # Search for spam similar to this post.
    similar = search_spam(text=text)

    # Add weight depending on number of post author has already made.
    # And the user score.
    authored = post.author.post_set.count()
    # Scales the values.

    # Get the weighted mean of a users activity score.
    #weighting_factor = (1 / (authored + post.author.profile.score))
    print("--")
    #print(weighting_factor, (authored + post.author.profile.score), weighting_factor)
    print("--")
    #weight = log(1 + (authored + post.author.profile.score) / (len([s.score for s in similar if s.is_spam]) + 5))

    scores = [s.score for s in similar if s.is_spam]

    # Take two local maximums and compute the n between them
    #nmaxs = len(scores)
    scores = sorted(scores, reverse=True)
    #scores = scores[:nmaxs]
    #scores = [s * weighting_factor for s in scores]
    #for m in range(nmaxs):
    #local_max = max(scores)
    #while local_max not in maxs and (nmaxs < len(maxs)):

    # Return the mean of the scores.
    if scores:
        mean = sum(scores) / len(scores)
    else:
        mean = 0

    #print(mean)
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


def test_classify(threshold=None, iterat=100):

    if threshold is None:
        threshold = settings.SPAM_THRESHOLD


    # Add posts to test spam index, then
    #test = add_file_to_index(fname=fname, delim=delim)

    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q())
    # Get the valid posts and shuffle.
    ham = Post.objects.valid_posts().order_by("?")

    ham = ham[:500]

    # Identify spam correctly.
    true_pos = 0

    # Identify valid post correctly.
    true_neg = 0

    # Missed to identify a spam.
    false_neg = 0

    # Mis-identified valid post as spam.
    false_pos = 0


    #for

    #
    # for spam in test:
    #     post_score = compute_score(text=spam)
    #
    #     return


    for post in chain([], ham):

        # Search for spam similar to this post.
        post_score = compute_score(text=post.content, post=post)

        # Positive spam

        if post_score >= threshold:
            if post.is_spam or post.is_deleted:
                #print(post_score, post.is_spam)
                #print(post.content)
                #print(post.is_spam or post.is_deleted, "TRUE")
                true_pos += 1
            else:
                #print(post_score, post.is_spam)
                #print(post.content)
                #print("FALSE ", post.is_spam or post.is_deleted)

                false_pos += 1
        # Negative spam
        else:
            if post.is_spam or post.is_deleted:
                #print("MISSED", post.is_spam or post.is_deleted)
                false_neg += 1
            else:
                true_neg += 1
            pass
        #print("-----")
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
    post_score = compute_score(text=post.text, post=post)

    # Update the spam score.
    Post.objects.filter(id=post.id).update(spam_score=post_score)

    # If the score exceeds threshold it gets quarantined.
    if post_score >= threshold:
        Post.objects.filter(id=post.id).update(spam=Post.SUSPECT)
        auth.log_action(log_text=f"Quarantined post={post.uid}; spam score={post_score}")
