import logging
import os
import time
from itertools import count, islice
from collections import defaultdict

# Postgres specific queries should go into separate module.
from django.conf import settings
from django.db.models import Q
from whoosh import writing, classify
from whoosh.analysis import StemmingAnalyzer
from whoosh.writing import AsyncWriter, BufferedWriter
from whoosh.searching import Results
import html2markdown
import bleach
from whoosh.qparser import MultifieldParser, OrGroup
from whoosh.analysis import STOP_WORDS
from whoosh.index import create_in, open_dir, exists_in
from whoosh.fields import ID, TEXT, KEYWORD, Schema, BOOLEAN, NUMERIC, DATETIME

from biostar.forum.models import Post

logger = logging.getLogger('engine')

# Stop words ignored where searching.
STOP = ['there', 'where', 'who', 'that'] + [w for w in STOP_WORDS]
STOP = set(STOP)


def timer_func():
    """
    Prints progress on inserting elements.
    """

    last = time.time()

    def elapsed(msg):
        nonlocal last
        now = time.time()
        sec = round(now - last, 1)
        last = now
        logger.debug(f"{msg} in {sec} seconds")

    def progress(index, step=500, total=0, msg=""):
        nonlocal last
        if index % step == 0:
            percent = int((index / total) * 100) if total >= index else index
            logger.debug(f"... {percent}% ({index} out of {total}). {step} {msg}")

    return elapsed, progress


def copy_hits(result, highlight=False):
    """
    Copy the items in results into a dict.
    """
    # Highlight title and content
    if highlight:
        title = result.highlights("title", minscore=0)
        title = title or result.get('title')
        content = result.highlights("content", top=5, minscore=0)
        content = content or result.get('content')
    else:
        title = result.get('title')
        content = result.get('content')

    bunched = dict(title=title,
                   content=content,
                   uid=result.get('uid'),
                   tags=result.get('tags'),
                   author=result.get('author'))
    return bunched


def index_exists(dirname=settings.INDEX_DIR, indexname=settings.INDEX_NAME):
    return exists_in(dirname=dirname, indexname=indexname)


def add_index(post, writer):

    # Ensure the content is stripped of any html.
    content = bleach.clean(post.content, styles=[], attributes={}, tags=[], strip=True)
    writer.update_document(title=post.title,
                           content=content,
                           tags=post.tag_val,
                           author=post.author.profile.name,
                           uid=post.uid)


def get_schema():
    analyzer = StemmingAnalyzer(stoplist=STOP)
    schema = Schema(title=TEXT(analyzer=analyzer, stored=True, sortable=True),
                    content=TEXT(analyzer=analyzer, stored=True, sortable=True),
                    tags=KEYWORD(commas=True, stored=True),
                    author=TEXT(stored=True),
                    uid=ID(stored=True))
    return schema

#
# def get_schema():
#    """
#    This is the issue!!!!!!
#    """
#     analyzer = StemmingAnalyzer(stoplist=STOP)
#     schema = Schema(title=TEXT(analyzer=analyzer, sortable=True),
#                     url=ID(stored=True),
#                     content_length=NUMERIC(sortable=True),
#                     thread_votecount=NUMERIC(stored=True, sortable=True),
#                     vote_count=NUMERIC(stored=True, sortable=True),
#                     content=TEXT(stored=True, analyzer=analyzer, sortable=True),
#                     tags=KEYWORD(stored=True, commas=True),
#                     is_toplevel=BOOLEAN(stored=True),
#                     author_is_moderator=BOOLEAN(stored=True),
#                     lastedit_user_is_moderator=BOOLEAN(stored=True),
#                     lastedit_user_is_suspended=BOOLEAN(stored=True),
#                     author_is_suspended=BOOLEAN(stored=True),
#                     lastedit_date=DATETIME(stored=True, sortable=True),
#                     creation_date=DATETIME(stored=True, sortable=True),
#                     rank=NUMERIC(stored=True, sortable=True),
#                     author=TEXT(stored=True),
#                     lastedit_user=TEXT(stored=True),
#                     lastedit_user_email=TEXT(stored=True),
#                     lastedit_user_score=NUMERIC(stored=True, sortable=True),
#                     lastedit_user_uid=ID(stored=True),
#                     lastedit_user_url=ID(stored=True),
#                     author_score=NUMERIC(stored=True, sortable=True),
#                     author_handle=TEXT(stored=True),
#                     author_email=TEXT(stored=True),
#                     author_uid=ID(stored=True),
#                     author_url=ID(stored=True),
#                     root_has_accepted=BOOLEAN(stored=True),
#                     reply_count=NUMERIC(stored=True, sortable=True),
#                     view_count=NUMERIC(stored=True, sortable=True),
#                     answer_count=NUMERIC(stored=True, sortable=True),
#                     uid=ID(stored=True),
#                     type=NUMERIC(stored=True, sortable=True),
#                     type_display=TEXT(stored=True))
#     return schema


def init_index(dirname=None, indexname=None, schema=None):
    # Initialize a new index or return an already existing one.

    ix_scheme = schema or get_schema()
    dirname = dirname or settings.INDEX_DIR
    indexname = indexname or settings.INDEX_NAME

    if exists_in(dirname=dirname, indexname=indexname):
        ix = open_dir(dirname=dirname, indexname=indexname)
    else:
        # Ensure index directory exists.
        os.makedirs(dirname, exist_ok=True)
        ix = create_in(dirname=dirname, schema=ix_scheme, indexname=indexname)

    return ix


def print_info(dirname=None, indexname=None):
    """
    Prints information on the index.
    """
    dirname = dirname or settings.INDEX_DIR
    indexname = indexname or settings.INDEX_NAME
    ix = init_index(dirname=dirname, indexname=indexname)

    counter = defaultdict(int)
    for index, fields in enumerate(ix.searcher().all_stored_fields()):
        key = fields['type_display']
        counter[key] += 1

    total = 0
    print('-' * 20)
    for key, value in counter.items():
        total += value
        print(f"{value}\t{key}")
    print('-' * 20)
    print(f"{total} total posts")


def index_posts(posts, ix=None, overwrite=False, add_func=add_index):
    """
    Create or update a search index of posts.
    """

    ix = ix or init_index()
    # The writer is asynchronous by default
    writer = AsyncWriter(ix)

    elapsed, progress = timer_func()
    total = posts.count()
    stream = islice(zip(count(1), posts), None)

    # Loop through posts and add to index
    for step, post in stream:
        progress(step, total=total, msg="posts indexed")
        add_func(post=post, writer=writer)

    # Commit to index
    if overwrite:
        logger.info("Overwriting the old index")
        writer.commit(mergetype=writing.CLEAR)
    else:
        logger.debug("Committing to index")
        writer.commit()

    elapsed(f"Committed {total} posts to index.")


def crawl(reindex=False, overwrite=False, limit=1000):
    """
    Crawl through posts in batches and add them to index.
    """

    if reindex:
        logger.info(f"Setting indexed field to false on all post.")
        Post.objects.filter(indexed=True).exclude(root=None).update(indexed=False)

    # Index a limited number of posts at one time.
    posts = Post.objects.valid_posts().exclude(Q(spam=Post.SPAM) | Q(indexed=False))[:limit]

    try:
        # Add post to search index.
        index_posts(posts=posts, overwrite=overwrite)
    except Exception as exc:
        logger.error(f'Error updating index: {exc}')
        Post.objects.filter(id__in=posts.values('id')).update(indexed=False)

    # Set the indexed field to true
    Post.objects.filter(id__in=posts.values('id')).update(indexed=True)

    return


def whoosh_search(query, limit=10, ix=None, fields=None, sortedby=[], **kwargs):
    """
    Query search index
    """

    fields = fields or ['tags', 'title', 'content', 'author']
    ix = ix or init_index()
    searcher = ix.searcher()

    # Splits the query into words and applies
    # and OR filter, eg. 'foo bar' == 'foo OR bar'
    orgroup = OrGroup

    parser = MultifieldParser(fieldnames=fields, schema=ix.schema, group=orgroup).parse(query)

    hits = searcher.search(parser,
                           sortedby=sortedby,
                           terms=True,
                           limit=limit)

    # Allow larger fragments
    hits.fragmenter.maxchars = 100
    hits.fragmenter.surround = 100

    return hits


def perform_search(query, fields=None, sortedby=[], limit=None):
    """
    Utility functions to search whoosh index, collect results and closes
    """

    limit = limit or settings.SEARCH_LIMIT

    hits = whoosh_search(query=query,
                         fields=fields,
                         sortedby=sortedby,
                         limit=limit)

    # Highlight the whoosh results.
    copier = lambda r: copy_hits(r, highlight=True)

    final = list(map(copier, hits))
    # Ensure searcher object gets closed.
    hits.searcher.close()

    return final


def more_like_this(uid, top=0, sortedby=[]):
    """
    Return posts in search index most similar to given post.
    """

    top = top or settings.SIMILAR_FEED_COUNT
    fields = ['uid']
    results = whoosh_search(query=uid, sortedby=sortedby, fields=fields)

    if len(results):
        results = results[0].more_like_this("content", top=top)
        # Copy hits to list and close searcher object.
        final = list(map(copy_hits, results))
        # Show unique posts
        final = {h['uid']: h for h in final}
        final = list(final.values())
    else:
        final = []

    # Ensure searcher object gets closed.
    results.searcher.close()

    return final


def remove_post(post, ix=None):
    """
    Remove spam from index
    """

    ix = ix or init_index()

    # Remove this post from index
    writer = AsyncWriter(ix)
    writer.delete_by_term('uid', text=post.uid)
    writer.commit()
    logger.debug(f"Removing uid={post.uid} from index")
    return
