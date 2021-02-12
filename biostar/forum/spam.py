import logging
import os
import shutil
import random
import time
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


def add_spam(uid):

    post = Post.objects.filter(uid=uid).first()
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


def search_spam(post, ix,):
    """
    Search spam index for posts similar to this one.
    Returns
    """
    writer = AsyncWriter(ix)
    add_post_to_index(post=post, writer=writer, is_spam=post.is_spam)
    writer.commit()

    # Search for this post in the spam index
    fields = ['uid']

    results = search.preform_whoosh_search(ix=ix, query=post.uid, fields=fields)

    # Preform more_like_this on this posts content
    similar_content = results[0].more_like_this('content', top=5)

    # Remove this post from the spam index after results are collected.
    writer = AsyncWriter(ix)
    writer.delete_by_term('uid', text=post.uid)
    writer.commit()

    # Get the results into a list and close the searcher object.
    similar_content = list(map(search.normalize_result, similar_content))

    results.searcher.close()

    return similar_content


def compute_score(post, ix=None):

    ix = ix or init_spam_index()
    N = 1
    weight = .7
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
    n = len(scores)

    scores = [s for s in scores][:N]

    # Return the mean of the scores.
    if scores:
        mean = abs(sum(scores) / len(scores[:N]))
    else:
        # Apply a weighted version of the threshold
        # when no similar posts are found.
        mean = settings.SPAM_THRESHOLD

    mean = ((mean * n) * weight) + bias

    return mean


def accuracy(tp, tn, fp, fn):

    return (tp + tn) / (tp + tn + fp + fn)


def specificity(tn, fp):
    return tn / (tn + fp)


def sensitivity(tp, fn):

    return tp / (tp + fn)


def sizer(lst, size):
    return size if size <= len(lst) else len(lst)


def false_positive_rate(fp, tn):
    return fp / (fp + tn)


def one_out_train(spam, ham, writer, size=100):

    spam = random.sample(spam, k=sizer(spam, size=size))
    ham = random.sample(ham, k=sizer(ham, size=size))

    # Remove one random item from the list
    one_out = random.choice([ham[0], spam[0]])
    ham, spam = ham[1:], spam[1:]

    # Get final queryset and only required fields
    posts = Post.objects.filter(id__in=chain(spam, ham)).exclude(id=one_out)
    posts = posts.only("content", "uid", "spam", "status", "title")

    # Write spam and ham posts to train index
    for post in posts:
        add_post_to_index(post=post, writer=writer, is_spam=post.is_spam)

    # Return the test post excluded from training.
    one_out = Post.objects.filter(id=one_out).first()
    return one_out


def report(nham, nspam, tn, tp, fn, fp):
    percent = lambda x: f"{x * 100:0.3f} %"
    acc = percent(accuracy(tp=tp, tn=tn, fp=fp, fn=fn))
    specif = percent(specificity(tn=tn, fp=fp))
    sens = percent(sensitivity(tp=tp, fn=fn))
    fps = percent(false_positive_rate(fp=fp, tn=tn))

    print(f"... \t{nspam}\tSPAM actual")
    print(f"... \t{tp}\tSPAM predicted\n\t\t---")
    print(f"... \t{nham}\tHAM actual")
    print(f"... \t{tn}\tHAM predicted")
    print("-"*10)
    print(f"... {acc}\tAccuracy\ttp + tn / (tp + tn + fp + fn) ")
    print(f"... {specif}\tSpecificity\ttn / (tn + fp)")
    print(f"... {sens}\tSensitivity\ttp / (tp + fn)")
    print(f"... {fps}\tFalse positive rate\tfp / (fp + tn)")
    return


def detail(post, post_score, is_spam=True, predict=True, verb=1):
    fp = (not is_spam) and predict
    fn = is_spam and (not predict)

    if verb and not (fn or fp):
        return
    print()
    if fn:
        print(f"-----\tFALSE NEGATIVE ( missed spam )\tuid={post.uid} score={post_score}. deleted={post.is_deleted}")

    elif fp:
        print(f"-----\tFALSE POSITIVE ( missed ham )\tuid={post.uid} score={post_score}. ")

    if verb > 1 and (fp or fn):
        print(post.content)
        print(">"*5)

    if verb >= 1 and (fp or fn):
        print("POSER SCORE", post_score)

    print()
    print("-" * 5)
    return


def test_classify(threshold=None, niter=100, limitmb=1024, size=100, verbosity=0):

    if threshold is None:
        threshold = settings.SPAM_THRESHOLD

    # Add posts to test spam index, then
    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED))

    # Get the valid posts and shuffle.
    ham = Post.objects.valid_posts(author__profile__score__lte=0, type__in=[Post.ANSWER, Post.COMMENT])

    # Get list of id's for both
    spam = list(spam.values_list("id", flat=True))
    ham = list(ham.values_list("id", flat=True))

    # tp = Identify spam correctly.
    # tn = Identify valid post correctly.
    # fn = Missed to identify a spam.
    # fp = Mis-identified valid post as spam.
    tp, tn, fn, fp = 0, 0, 0, 0
    seen_ham, seen_spam = 0, 0
    elapsed, progress = util.timer_func()

    for i in range(niter):
        # Remove previous index
        if os.path.exists(TRAIN_DIR):
            shutil.rmtree(TRAIN_DIR)

        ix = search.init_index(dirname=TRAIN_DIR,
                               indexname=f"train_{util.get_uuid(8)}_{settings.SPAM_INDEX_NAME}",
                               schema=spam_schema())
        writer = BufferedWriter(ix, limit=int((niter/2) + 1), writerargs=dict(limitmb=limitmb, multisegment=True))

        index_writer(writer=writer, title="Placeholder",
                     content_length=0, is_spam=True,
                     content='CONTENT', uid=STARTER_UID)

        # Take one spam post out of training set.
        one_out = one_out_train(spam=spam, writer=writer, size=size, ham=ham)
        writer.commit()
        writer.close()
        post_score = compute_score(post=one_out, ix=ix)

        predicted_spam = post_score >= threshold
        is_spam = one_out.is_spam or one_out.is_deleted
        is_ham = not is_spam

        seen_spam += 1 if is_spam else 0
        seen_ham += 1 if is_ham else 0

        detail(is_spam=is_spam, predict=predicted_spam, post=one_out,
               verb=verbosity, post_score=post_score)

        if predicted_spam:
            tp += 1 if is_spam else 0
            fp += 1 if is_ham else 0

        else:
            fn += 1 if is_spam else 0
            tn += 1 if is_ham else 0

        progress(i, step=5, msg=f"iterations. tp={tp} fp={fp} tn={tn} fn={fn}")

    train_spam = sizer(spam, size)
    train_ham = sizer(ham, size)
    print(f"... {train_ham + train_spam}\tSize of index ( per iteration )")
    print(f"... \t{train_spam}\tSPAM")
    print(f"... \t{train_ham}\tHAM")
    print(f"\n... {niter}\tNumber of iterations")
    report(nham=seen_ham, nspam=seen_spam, tn=tn, tp=tp, fp=fp, fn=fn)

    return


def score(uid, threshold=None):
    """
    """

    if not settings.CLASSIFY_SPAM:
        return

    if threshold is None:
        threshold = settings.SPAM_THRESHOLD

    post = Post.objects.filter(uid=uid).first()

    # User's with high enough score automatically given green light.
    # if not post.author.profile.low_rep:
    #    return

    # Search for spam similar to this post.
    post_score = compute_score(post=post)

    # Update the spam score.
    Post.objects.filter(id=post.id).update(spam_score=post_score)

    # If the score exceeds threshold it gets labeled spam.
    if post_score >= threshold:
        Post.objects.filter(id=post.id).update(spam=Post.SPAM)
        auth.db_logger(text=f" post={post.uid}; spam score={post_score}")
