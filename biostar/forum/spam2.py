import logging
import random
from itertools import chain
from django.conf import settings
from django.db.models import Q
from elasticsearch_dsl import Search, Index
from elasticsearch_dsl import Q as elastic_Q
from elasticsearch_dsl.query import MoreLikeThis, Query
from elasticsearch_dsl.connections import connections
from elasticsearch_dsl import FacetedSearch, TermsFacet, DateHistogramFacet
from elasticsearch import Elasticsearch
from elasticsearch.helpers import bulk
from biostar.forum.models import Post
from biostar.forum import util, search
from biostar.forum.documents import SpamDocument, TrainSpam

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


def to_dict(posts):

    data = ({
        "_index": "train_spam",
        "_id": p.id,
        "_title": p.title,
        "_content": p.content,
        "_is_spam": p.is_spam or p.is_deleted
    }
    for p in posts.iterator()
    )
    return data


def report(nham, nspam, tn, tp, fn, fp):
    percent = lambda x: f"{x * 100:0.3f} %"

    acc = percent(accuracy(tp=tp, tn=tn, fp=fp, fn=fn))
    specif = percent(specificity(tn=tn, fp=fp))
    sens = percent(sensitivity(tp=tp, fn=fn))
    missed = percent(miss_rate(fn=fn, tp=tp))
    false_discovery = percent(false_discovery_rate(fp=fp, tp=tp))

    print(f"Total data in each iteration (spam + ham ) : {nspam + nham}")
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


def search_spam(post, indexname=None, client=None):

    s = Search(using=client, index=indexname, doc_type=TrainSpam)

    text = [post.content.replace(f" {stop} ", "") for stop in search.STOP]

    res = s.suggest("suggest_1", text=text, phrase=dict(field='content',
                                                        size=1,
                                                        confidence=0,
                                                        real_word_error_likelihood=1))
    res = res.execute()

    res8 = [x.text for x in res.suggest.suggest_1]
    #s.query()

    print(len(res8), type(res.suggest), res.hits, post.content, post.is_spam, res.to_dict())
    1/0
    return s


def score(post, indexname=None, client=None):

    client = client or Elasticsearch()

    similar = search_spam(client=client, post=post, indexname=indexname)

    print([x for x in similar], post.is_spam)
    #similar = search(post=post, more_like_this=True)
    1/0

    return 0


def bulk_index(qs, doc=None):

    doc = doc or SpamDocument
    client = Elasticsearch()

    if client.indices.exists(index=doc.Index.name):
        client.indices.delete(index=doc.Index.name, ignore=[400, 404])
        logger.info("Created indexed.")
        doc.init()

    print(len([d for d in to_dict(qs)]))

    status = bulk(client=client, actions=to_dict(qs))

    success, fail = 1, 0
    client.indices.refresh(doc.Index.name)
    #print(client.indices.refresh(doc.Index.name))
    print(client.cat.count(doc.Index.name, params={"format": "json"}))
    if status == fail:
        logger.error(f"Error building index.")


def train(spamids, hamids, limit=None):
    limit = limit or 500

    random.shuffle(spamids)
    random.shuffle(hamids)

    one_out = [spamids[0], hamids[0]]
    one_out = random.choice(one_out)

    spamids, hamids = spamids[1:limit], hamids[1:limit]

    qs = Post.objects.filter(id__in=chain(hamids, spamids))
    one_out = Post.objects.filter(id=one_out).first()

    bulk_index(qs=qs, doc=TrainSpam)

    return one_out


def test(threshold=None, niter=100, limit=None):

    threshold = threshold or settings.SPAM_THRESHOLD

    indexname = settings.TRAIN_SPAM_INDEX

    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED))
    ham = Post.objects.valid_posts()

    spamids = list(spam.values_list('id', flat=True))
    hamids = list(ham.values_list('id', flat=True))

    _, progress = util.timer_func()

    # tp = Identify spam correctly.
    # tn = Identify valid post correctly.
    # fn = Missed to identify a spam.
    # fp = Mis-identified valid post as spam.
    tp, tn, fn, fp = 0, 0, 0, 0

    for n in range(niter):
        one_out = train(spamids=spamids, hamids=hamids, limit=limit)
        is_spam = one_out.is_spam or one_out.is_deleted
        spam_score = score(post=one_out, indexname=indexname)

        if spam_score >= threshold:
            tp += 1 if is_spam else 0
            fp += 1 if (not is_spam) else 0
        else:
            fn += 1 if is_spam else 0
            tn += 1 if (not is_spam) else 0

        progress(n, step=10, msg=f"iterations; tp={tp} fp={fp} tn={tn} fn={fn}")

    nham = len(spamids[:limit])
    nspam = len(hamids[:limit])

    report(nham=nham, nspam=nspam, tn=tn, fp=fp, tp=tp, fn=fn)


