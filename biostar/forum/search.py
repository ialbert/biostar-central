import logging
import os
import time
from itertools import count, islice
from collections import defaultdict

# Postgres specific queries should go into separate module.
from django.conf import settings
from django.db.models import Q
from whoosh import writing, classify
from whoosh.analysis import StemmingAnalyzer, StopFilter
from whoosh.writing import AsyncWriter, BufferedWriter
from whoosh.searching import Results, ResultsPage
import html2markdown
import bleach
from whoosh.qparser import MultifieldParser, OrGroup
from whoosh.analysis import STOP_WORDS
from whoosh.index import create_in, open_dir, exists_in
from whoosh.fields import ID, TEXT, KEYWORD, Schema, BOOLEAN, NUMERIC, DATETIME

from biostar.utils.helpers import htmltomarkdown
from biostar.forum.models import Post

logger = logging.getLogger('engine')

# Stop words ignored where searching.
STOP = ['there', 'where', 'who', 'that', 'to', 'do', 'my', 'only', 'but', 'about',
        'our', 'able', 'how', 'am', 'so', 'want']

STOP += [w for w in STOP_WORDS]
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
        title = result.highlights("title", minscore=0) or result.get('title')
        content = result.highlights("content", top=5, minscore=0) or result.get('content')
        tags = result.highlights("tags", minscore=0) or result.get('tags')
    else:
        title = result.get('title')
        content = result.get('content')
        tags = result.get('tags')

    bunched = dict(title=title,
                   content=content,
                   uid=result.get('uid'),
                   tags=tags,
                   author=result.get('author'),
                   lastedit_date=result.get('lastedit_date'))
    return bunched


def index_exists(dirname=settings.INDEX_DIR, indexname=settings.INDEX_NAME):
    return exists_in(dirname=dirname, indexname=indexname)


def add_index(post, writer):
    # Ensure the content is stripped of any html.
    content = htmltomarkdown(post.content)

    writer.update_document(title=post.title,
                           content=content,
                           tags=post.tag_val,
                           author=post.author.profile.name,
                           uid=post.uid,
                           lastedit_date=post.lastedit_date)


def get_schema():
    analyzer = StemmingAnalyzer(stoplist=STOP) | StopFilter(stoplist=STOP)
    schema = Schema(title=TEXT(analyzer=analyzer, stored=True, sortable=True),
                    content=TEXT(analyzer=analyzer, stored=True, sortable=True),
                    tags=KEYWORD(commas=True, stored=True),
                    author=TEXT(stored=True),
                    uid=ID(unique=True, stored=True),
                    lastedit_date=DATETIME(sortable=True, stored=True))
    return schema


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


def whoosh_search(query, limit=10, page=1, ix=None, fields=None, reverse=False, sortedby=[], **kwargs):
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

    hits = searcher.search_page(parser,pagenum=page, pagelen=limit, reverse=reverse, sortedby=sortedby, **kwargs)
    hits.results.fragmenter.maxchars = 100
    hits.results.fragmenter.surround = 100

    return hits


def perform_search(query, page=1, fields=None, reverse=False, sortedby=[], limit=None):
    """
    Utility functions to search whoosh index, collect results and closes
    """

    limit = limit or settings.SEARCH_LIMIT

    indexed = whoosh_search(query=query, fields=fields, page=page, reverse=reverse, sortedby=sortedby, limit=limit)

    # Highlight the whoosh results.
    copier = lambda r: copy_hits(r, highlight=True)

    final = list(map(copier, indexed))

    indexed.results.searcher.close()

    return final, indexed


def more_like_this(uid, top=0, sortedby=[]):
    """
    Return posts in search index most similar to given post.
    """

    top = top or settings.SIMILAR_FEED_COUNT
    fields = ['uid']
    found = whoosh_search(query=uid, sortedby=sortedby, fields=fields)

    if len(found):
        hits = found[0].more_like_this("content", top=top)
        # Copy hits to list and close searcher object.
        final = list(map(copy_hits, hits))
    else:
        final = []

    found.results.searcher.close()

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
