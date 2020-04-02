import logging
from itertools import groupby, islice, count, chain
from django.conf import settings
from django.db.models import Q
from whoosh.writing import AsyncWriter
from whoosh.analysis import StemmingAnalyzer
from whoosh.fields import ID, TEXT, KEYWORD, Schema, NUMERIC
from biostar.forum.models import Post
from biostar.forum import search, auth, util

logger = logging.getLogger("engine")


def spam_schema():
    analyzer = StemmingAnalyzer()
    schema = Schema(title=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    content=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    content_length=NUMERIC(stored=True, sortable=True),
                    uid=ID(stored=True),
                    is_spam=True)
    return schema


def index_writer(writer, **kwargs):
    writer.update_document(title=kwargs.get("title"),
                           content_length=kwargs.get("content_length"),
                           content=kwargs.get("content"),
                           uid=kwargs.get("uid"),
                           is_spam=kwargs.get("is_spam"))


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
                     content_length=len(text),
                     content=text, uid= util.get_uuid(16))

    writer.commit()
    elapsed("Committed to index.")


def add_post_to_index(post, writer):
    """
    Insert post to index.
    """
    index_writer(writer=writer, title=post.title,
                 content_length=len(post.content),
                 content=post.content,
                 uid=post.uid)


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


def compute_score(post):
    # Search for spam similar to this post.
    similar = search_spam(post=post)

    # Get the similarity score for the returned spam
    scores = [s.score for s in similar if s.is_spam]

    # Return the mean of the scores.
    if scores:
        mean = sum(scores) / len(scores)
    else:
        mean = 0

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
    scount = spam.count()

    # Get the valid posts.
    ham = Post.objects.filter(Q(spam=Post.NOT_SPAM) | Q(status=Post.DEFAULT))[:scount]

    # Identify spam correctly.
    true_pos = 0

    # Identify valid post correctly.
    true_neg = 0

    # Missed to identify a spam.
    false_neg = 0

    # Mis-identified valid post as spam.
    false_pos = 0

    for post in chain(spam, ham):

        # Search for spam similar to this post.
        post_score = compute_score(post=post)

        #if not post.author.profile.low_rep:
        #    continue

        # Positive spam
        if post_score >= threshold:
            if post.is_spam or post.is_deleted:
                print(post_score, post.is_spam)
                print(post.content)
                true_pos += 1
            else:
                false_pos += 1


        # Negative spam
        else:
            if post.is_spam or post.is_deleted:
                false_neg += 1
            else:
                true_neg += 1
            pass

    acc = accuracy(tp=true_pos, tn=true_neg, fp=false_pos, fn=false_neg)
    specif = specificity(tn=true_neg, fp=false_pos)
    sens = sensitivity(tp=true_neg, fn=false_neg)
    missed = miss_rate(fn=false_neg, tp=true_pos)
    false_discovery = false_discovery_rate(fp=false_pos, tp=true_pos)

    print(f"Number of ham: {ham.count()}")
    print(f"Number of spam: {spam.count()}")
    print(f"Accuracy: {acc} {true_pos}")
    print(f"Specificity: {specif}")
    print(f"Sensitivity: {sens}")
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
