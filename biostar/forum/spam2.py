import logging
import random
from itertools import chain
from django.conf import settings
from django.db.models import Q
from elasticsearch_dsl import Search
from elasticsearch_dsl import Q as elastic_Q
from elasticsearch import Elasticsearch
from elasticsearch.helpers import bulk
from biostar.forum.models import Post
from biostar.forum import util
from biostar.forum.documents import SpamDocument

SPAM_LIMIT = 2500
HAM_LIMIT = 10000

logger = logging.getLogger("engine")


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


def report(nham, nspam, tn, tp, fn, fp):
    percent = lambda x: f"{x * 100:0.3f} %"

    acc = percent(accuracy(tp=tp, tn=tn, fp=fp, fn=fn))
    specif = percent(specificity(tn=tn, fp=fp))
    sens = percent(sensitivity(tp=tp, fn=fn))
    missed = percent(percent(miss_rate(fn=fn, tp=tp)))
    false_discovery = false_discovery_rate(fp=fp, tp=tp)

    print(f"Actual ham: {nham}")
    print(f"Predicted ham: {tn}")
    print(f"Actual spam: {nspam}")
    print(f"Predicted spam: {tp}")

    print(f"... {acc}\tAccuracy\ttp + tn / (tp + tn + fp + fn) ")
    print(f"... {specif}\tSpecificity\ttn / (tn + fp)")
    print(f"... {sens}\tSensitivity\ttp / (tp + fn)")
    print(f"... {missed}\tMissed rate\tfn / (fn + tp) ")
    print(f"... {false_discovery}\tFalse discovery rate\t{false_discovery}")
    return


def search(doc, post, more_like_this=True):
    SpamDocument(indexname=settings.SPAM_INDEX_NAME).init()

    client = Elasticsearch()
    q = Q("bool", should=[elastic_Q("match", firstname=firstname),
                          elastic_Q("match", gender=gender)], minimum_should_match=1)
    s = Search(using=client, index="bank").query(q)[0:20]
    response = s.execute()
    #spam = doc()

    return


def score(post, indexname=None):

    #if indexname is None:

    similar = search(doc=SpamDocument, post=post, more_like_this=True)
    1/0

    return 0


def bulk_index(qs):

    es = Elasticsearch()
    elapsed, _ = util.timer_func()
    status = bulk(client=es, actions=(obj.spam_indexing() for obj in qs.iterator()))
    count = qs.count()
    success, fail = 1, 0

    if status == fail:
        logger.error(f"Error building index.")

    elapsed(f"Finished adding with {count} objects")


def train(indexname, spamids, hamids, limit=500):

    es = Elasticsearch()
    if es.indices.exists(index=indexname):
        es.indices.delete(index=indexname, ignore=[400, 404])
        logger.info("Removed indexed.")

    random.shuffle(spamids)
    random.shuffle(hamids)

    one_out = [spamids.pop(0), hamids.pop(0)]
    one_out = random.choice(one_out)

    hamids, spamids = hamids[:limit], spamids[:limit]

    qs = Post.objects.filter(id__in=chain(hamids, spamids))
    one_out = Post.objects.filter(id=one_out).first()

    bulk_index(qs=qs)

    return one_out


def test(indexname, threshold=None, niter=100, limit=500):

    threshold = threshold or settings.SPAM_THRESHOLD
    SpamDocument(indexname=indexname).init()
    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED))
    ham = Post.objects.valid_posts()

    spamids = spam.values_list('id', flat=True)
    hamids = ham.values_list('id', flat=True)
    _, progress = util.timer_func()

    # tp = Identify spam correctly.
    # tn = Identify valid post correctly.
    # fn = Missed to identify a spam.
    # fp = Mis-identified valid post as spam.
    tp, tn, fn, fp = 0, 0, 0, 0

    for n in range(niter):

        one_out = train(indexname=indexname, spamids=spamids, hamids=hamids, limit=limit)
        is_spam = one_out.is_spam or one_out.is_deleted
        spam_score = score(post=one_out, indexname=indexname)

        if spam_score >= threshold:
            tp += 1 if is_spam else 0
            fp += 1 if (not is_spam) else 0
        else:
            fn += 1 if is_spam else 0
            tn += 1 if (not is_spam) else 0

        progress(n, step=100, msg=f"iterations; tp={tp} fp={fp} tn={tn} fn={fn}")

    nham = len(spamids[:limit])
    nspam = len(hamids[:limit])

    report(nham=nham, nspam=nspam, tn=tn, fp=fp, tp=tp, fn=fn)


# def build_index():
#
#     article = SpamDocument(meta={'id': 42}, title='Hello world!', content="Test", is_spam=True)
#     article.save()
#
#     article = SpamDocument.get(id=42)
#     print(article.is_spam)
#     return

