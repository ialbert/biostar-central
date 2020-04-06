import logging
import os
import shutil
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


def add_post_to_index(post, writer, is_spam=False):
    """
    Insert post to index.
    """
    index_writer(writer=writer, title=post.title,
                 content_length=len(post.content),
                 content=post.content, is_spam=is_spam,
                 uid=f"{post.uid}")


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


def search_spam(post, ix, include_title=False):
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
    if include_title:
        similar_title = results[0].more_like_this('title', text=post.title, top=5000)

    # Preform more_like_this on this posts content
    similar_content = results[0].more_like_this('content', top=500)
    similar_title = []

    # Remove this post from the spam index after results are collected.
    writer = AsyncWriter(ix)
    writer.delete_by_term('uid', text=post.uid)
    writer.commit()

    # Get the results into a list and close the searcher object.
    similar_content = list(map(search.normalize_result, similar_content))

    #results.searcher.close()

    return similar_content, similar_title


def compute_score(post, ix=None):

    ix = ix or init_spam_index()

    # Search for spam similar to this post.
    similar_content, similar_title = search_spam(post=post, ix=ix, include_title=True)

    # Add weight depending on number of post author has already made.
    # And the user score.
    #authored = post.author.post_set.count()
    # Scales the values.

    # Get the weighted mean of a users activity score.

    weighting_factor = 0.6

    authored = post.author.post_set.exclude(id=post.id).count() * 5

    boost = (1 - weighting_factor) + 2 if (post.is_comment or post.is_job) else 0

    weighting_factor += boost

    weight = (1 / ((post.author.profile.score * 10 + authored) + weighting_factor))

    scores = [s.score for s in similar_content if s.is_spam]
    #title_scores = [h.score for h in similar_title if h.get("is_spam")]

    # Take two local maximums and compute the n between them
    N = 1
    scores = sorted(scores, reverse=True)
    scores = [s * weight for s in scores][:N]

    # Return the mean of the scores.
    if scores:
        mean = sum(scores) / len(scores)
    else:
        mean = weight - log(2) ** 3

    print(scores, mean, post.title, weight, post.author.profile.score)
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

    ham = random.sample(ham, k=sizer(ham, size=size))
    spam = random.sample(spam, k=sizer(spam, size=size))

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

    fp = not is_spam and predict
    fn = is_spam and not predict

    if verb and not (fn or fp):
        return

    if fn:
        print(f"-----\tFALSE NEGATIVE ( missed spam )\tuid={post.uid} score={post_score}. deleted={post.is_deleted}")

    elif fp:
        print(f"++++++\tFALSE POSITIVE ( missed ham )\tuid={post.uid} score={post_score}. ")

    if verb > 1 and (fp or fn):
        print(post.content)
        print(">"*5)
        print(post.author)

    print("-" * 5)
    return


def test_classify(threshold=None, niter=100, size=100, verbosity=0):

    if threshold is None:
        threshold = settings.SPAM_THRESHOLD

    # Add posts to test spam index, then
    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED))

    # Get the valid posts and shuffle.
    ham = Post.objects.valid_posts()

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
        shutil.rmtree(TRAIN_DIR)

        ix = search.init_index(dirname=TRAIN_DIR, indexname=f"train_{util.get_uuid(8)}_{settings.SPAM_INDEX_NAME}",
                               schema=spam_schema())
        writer = BufferedWriter(ix, limit=niter, writerargs=dict(limitmb=1026, multisegment=True))

        index_writer(writer=writer, title="Placeholder",
                     content_length=0, is_spam=True,
                     content='CONTENT', uid=STARTER_UID)

        # Take one spam post out of training set.
        one_out = one_out_train(ham=ham, spam=spam, writer=writer, size=size)
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
