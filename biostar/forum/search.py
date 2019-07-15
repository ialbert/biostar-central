import logging
import time
import os
from itertools import count, islice

from django.conf import settings
from django.template import loader

from whoosh import writing
from whoosh.sorting import FieldFacet, ScoreFacet
from whoosh.writing import AsyncWriter
from whoosh.qparser import MultifieldParser
from whoosh.analysis import SpaceSeparatedTokenizer, StopFilter, STOP_WORDS
from whoosh.index import create_in, open_dir, exists_in
from whoosh.fields import ID, NGRAM, TEXT, KEYWORD, Schema, BOOLEAN, NUMERIC, NGRAMWORDS

logger = logging.getLogger('engine')


# Stop words ignored where searching.
STOP = ['there', 'where', 'who'] + [w for w in STOP_WORDS]
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
        print(f"{msg} in {sec} seconds")

    def progress(index, step=500, msg=""):
        nonlocal last
        if index % step == 0:
            elapsed(f"... {index} {msg}")

    return elapsed, progress


def index_exists():
    return exists_in(dirname=settings.INDEX_DIR, indexname=settings.INDEX_NAME)


def add_index(post, writer):
    writer.update_document(title=post.title, url=post.get_absolute_url(),
                           type_display=post.get_type_display(),
                           type=post.type,
                           lastedit_date=post.lastedit_date.timestamp(),
                           content=post.content, tags=post.tag_val,
                           is_toplevel=post.is_toplevel,
                           rank=post.rank, uid=post.uid,
                           author_handle=post.author.username,
                           author=post.author.profile.name,
                           author_score=post.author.profile.score,
                           author_uid=post.author.profile.uid,
                           author_url=post.author.profile.get_absolute_url())


def get_schema():
    tokenizer = SpaceSeparatedTokenizer() | StopFilter(stoplist=STOP)
    schema = Schema(title=NGRAMWORDS(stored=True, tokenizer=tokenizer, sortable=True),
                    url=ID(stored=True),
                    content=NGRAMWORDS(stored=True, tokenizer=tokenizer),
                    tags=KEYWORD(stored=True),
                    is_toplevel=BOOLEAN(stored=True),
                    lastedit_date=NUMERIC(stored=True, sortable=True),
                    rank=NUMERIC(stored=True, sortable=True),
                    author=TEXT(stored=True),
                    author_score=NUMERIC(stored=True, sortable=True),
                    author_handle=TEXT(stored=True),
                    author_uid=ID(stored=True),
                    author_url=ID(stored=True),
                    uid=ID(stored=True),
                    type=NUMERIC(stored=True, sortable=True),
                    type_display=TEXT(stored=True),)
    return schema


def init_index():
    # Initialize a new index or return an already existing one.

    if index_exists():
        ix = open_dir(dirname=settings.INDEX_DIR, indexname=settings.INDEX_NAME)
    else:
        # Ensure index directory exists.
        os.makedirs(settings.INDEX_DIR, exist_ok=True)
        ix = create_in(dirname=settings.INDEX_DIR, schema=get_schema(), indexname=settings.INDEX_NAME)

    return ix


def index_posts(posts, reindex=False):
    """
    Create or update a search index of posts.
    """

    ix = init_index()
    # The writer is asynchronous by default
    writer = AsyncWriter(ix)

    elapsed, progress = timer_func()
    stream = islice(zip(count(1), posts), None)

    # Loop through posts and add to index
    for step, post in stream:
        progress(step, msg="posts indexed")

        add_index(post=post, writer=writer)

    # Commit to index
    if reindex:
        # Re-index posts from scratch when committing.
        writer.commit(mergetype=writing.CLEAR)
    else:
        writer.commit()

    elapsed(f"""Indexed {len(posts)} posts: 
            dir={settings.INDEX_DIR} name={settings.INDEX_NAME}.""")


def query(q='', fields=['content'], **kwargs):
    """
    Query the indexed, looking for a match in the specified fields.
    Results a tuple of results and an open searcher object.
    """

    # Do not preform any queries if the index does not exist.
    if not index_exists():
        return []

    ix = init_index()
    searcher = ix.searcher()

    profile_score = FieldFacet("author_score", reverse=True)
    post_type = FieldFacet("type")
    rank = FieldFacet("rank", reverse=True)
    default = ScoreFacet()

    # Sort by: toplevel, match score, author reputation, post rank.
    sort_by = [post_type,  profile_score, rank, default]

    parser = MultifieldParser(fieldnames=fields, schema=ix.schema).parse(q)
    results = searcher.search(parser, sortedby=sort_by, limit=settings.SEARCH_LIMIT, terms=True, **kwargs)
    # Allow larger fragments
    results.fragmenter.maxchars = 100
    #results.fragmenter.charlimit = None
    # Show more context before and after
    results.fragmenter.surround = 100

    return results
