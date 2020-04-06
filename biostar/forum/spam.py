import logging
import os
import random
from math import log
from itertools import groupby, islice, count, chain
from django.conf import settings
from django.db.models import Q
from whoosh.writing import AsyncWriter, BufferedWriter
from whoosh import classify
from whoosh.analysis import StemmingAnalyzer
from whoosh.fields import ID, TEXT, KEYWORD, Schema, NUMERIC, BOOLEAN
from biostar.forum.models import Post
from biostar.forum import search, auth, util
from biostar.forum.documents import SpamDocument

logger = logging.getLogger("engine")

HAM_LIMIT = 5000


STARTER_UID = 'placeholder'


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


def add_post_to_index(post, writer, is_spam=False):
    """
    Insert post to index.
    """
    index_writer(writer=writer, title=post.title,
                 content_length=len(post.content),
                 content=post.content, is_spam=is_spam,
                 uid=post.uid)


def build_spam_index(overwrite=False):
    # Get all un-indexed spam posts.
    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED), indexed=False)

    # Initialize the spam index
    ix = bootstrap_index()

    # Index spam posts.
    search.index_posts(posts=spam, ix=ix, overwrite=overwrite, add_func=add_post_to_index)

    spam.update(indexed=True)
    logger.info("Built spam index.")
    return ix


def search_spam(post, ix):
    """
    Search spam index for posts similar to this one.
    Returns
    """
    writer = AsyncWriter(ix)
    add_post_to_index(post=post, writer=writer, is_spam=post.is_spam)
    writer.commit()

    # Search for this post in the spam index
    fields = ['uid']
    writer.searcher()
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


def compute_score(post, ix=None):

    ix = ix or init_spam_index()

    # Search for spam similar to this post.
    similar = search_spam(post=post, ix=ix)

    # Add weight depending on number of post author has already made.
    # And the user score.
    #authored = post.author.post_set.count()
    # Scales the values.

    # Get the weighted mean of a users activity score.
    weighting_factor = (1 / (post.author.profile.score + 1))
    # print("--")
    # print(weighting_factor, (authored + post.author.profile.score), weighting_factor)
    # print("--")

    scores = [s.score for s in similar if s.is_spam]

    # Take two local maximums and compute the n between them

    scores = sorted(scores, reverse=True)

    scores = [s * weighting_factor for s in scores]

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


def one_out_train(spam, ham, writer, size=100):

    # Build separate index then clean the index

    spam = spam[:size].values_list("id", flat=True)
    ham = ham[:size].values_list("id", flat=True)
    out_spam = random.choice(spam)
    out_ham = random.choice(ham)

    one_out = random.choice([out_ham, out_spam])
    posts = Post.objects.filter(id__in=chain(spam, ham)).exclude(id=one_out)
    posts = posts.only("content", "uid", "spam", "status", "title")

    # Write spam and hap posts to train index
    for post in posts:
        add_post_to_index(post=post, writer=writer, is_spam=post.is_spam)

    one_out = Post.objects.filter(id=one_out).first()
    # Return the test set excluded from training.
    return one_out


def report(tested_ham, tested_spam, tn, tp, fn, fp):
    acc = accuracy(tp=tp, tn=tn, fp=fp, fn=fn)
    specif = specificity(tn=tn, fp=fp)
    sens = sensitivity(tp=tp, fn=fn)
    missed = miss_rate(fn=fn, tp=tp)
    false_discovery = false_discovery_rate(fp=fp, tp=tp)

    print(f"Number of ham: {tested_ham}")
    print(f"Predicted ham: {tn}")
    print(f"Number of spam: {tested_spam}")
    print(f"Predicted spam: {tp}")

    print(f"Accuracy: {acc}")
    print(f"Specificity : {specif}")
    print(f"Sensitivity : {sens}")
    print(f"Missed rate: {missed}")
    print(f"False discovery rate:{false_discovery}")
    return


def test_classify(threshold=None, niter=100):

    if threshold is None:
        threshold = settings.SPAM_THRESHOLD

    # Add posts to test spam index, then
    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED))

    # Get the valid posts and shuffle.
    ham = Post.objects.valid_posts()

    # tp = Identify spam correctly.
    # tn = Identify valid post correctly.
    # fn = Missed to identify a spam.
    # fp = Mis-identified valid post as spam.
    tp, tn, fn, fp = 0, 0, 0, 0
    size = 90
    elapsed, progress = util.timer_func()
    tested_ham, tested_spam = 0, 0
    dirname = os.path.dirname(settings.SPAM_INDEX_DIR)
    basename = f"train_{os.path.basename(settings.SPAM_INDEX_DIR)}"
    train_dir = os.path.join(dirname, basename)
    import shutil

    for i in range(niter):
        shutil.rmtree(train_dir)
        ix = search.init_index(dirname=train_dir, indexname=f"train_{util.get_uuid(8)}_{settings.SPAM_INDEX_NAME}",
                               schema=spam_schema())
        writer = BufferedWriter(ix, limit=niter, writerargs=dict(limitmb=512, multisegment=True))

        index_writer(writer=writer, title="Placeholder",
                     content_length=0, is_spam=True,
                     content='CONTENT', uid=STARTER_UID)

        # Take one spam post out of training set.
        one_out = one_out_train(ham=ham, spam=spam, writer=writer, size=size)
        writer.commit()
        writer.close()

        post_score = compute_score(post=one_out, ix=ix)

        is_spam = one_out.is_spam or one_out.is_deleted
        is_ham = not is_spam

        tested_spam += 1 if is_spam else 0
        tested_ham += 1 if is_ham else 0

        if post_score >= threshold:
            tp += 1 if is_spam else 0
            fp += 1 if is_ham else 0
        else:
            if is_spam:
                print("FALSE", one_out.content,  one_out.is_spam, one_out.author.id,one_out.author)

            fn += 1 if is_spam else 0
            tn += 1 if is_ham else 0

        progress(i, step=5, msg=f"iterations. tp={tp} fp={fp} tn={tn} fn={fn}")

    elapsed(f"Results gathered over {niter} iterations.")
    report(tested_ham=tested_ham, tested_spam=tested_spam, tn=tn, tp=tp, fp=tp, fn=fn)

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
    post_score = compute_score(post=post)

    # Update the spam score.
    Post.objects.filter(id=post.id).update(spam_score=post_score)

    # If the score exceeds threshold it gets quarantined.
    if post_score >= threshold:
        Post.objects.filter(id=post.id).update(spam=Post.SUSPECT)
        auth.log_action(log_text=f"Quarantined post={post.uid}; spam score={post_score}")
