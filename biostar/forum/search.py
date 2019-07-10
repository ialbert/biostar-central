import logging
import time
import os
from itertools import count, islice

from django.conf import settings

from whoosh import writing
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

    def progress(index, step=100, msg=""):
        nonlocal last
        if index % step == 0:
            elapsed(f"... {index} {msg}")

    return elapsed, progress


def index_exists():
    return exists_in(dirname=settings.INDEX_DIR, indexname=settings.INDEX_NAME)


def add_index(post, writer):
    title = '' if not post.is_toplevel else post.title
    writer.add_document(title=title, url=post.get_absolute_url(),
                        type=post.get_type_display(),
                        lastedit_date=post.lastedit_date.timestamp(),
                        content=post.content, tags=post.tag_val,
                        is_toplevel=post.is_toplevel,
                        rank=post.rank, uid=post.uid,
                        author=post.author.profile.name,
                        author_uid=post.author.profile.uid,
                        author_url=post.author.profile.get_absolute_url())


def get_schema():
    tokenizer = SpaceSeparatedTokenizer() | StopFilter(stoplist=STOP)

    schema = Schema(title=NGRAMWORDS(stored=True, tokenizer=tokenizer),
                    url=ID(stored=True),
                    content=NGRAMWORDS(stored=True, tokenizer=tokenizer),
                    tags=KEYWORD(stored=True),
                    is_toplevel=BOOLEAN(stored=True),
                    lastedit_date=NUMERIC(stored=True),
                    author_uid=ID(stored=True),
                    rank=NUMERIC(stored=True, sortable=True),
                    author=TEXT(stored=True),
                    author_url=ID(stored=True),
                    uid=ID(stored=True),
                    type=TEXT(stored=True))
    return schema


def init_index():
    # Initialize index or return already existing one.
    if index_exists():
        ix = open_dir(dirname=settings.INDEX_DIR, indexname=settings.INDEX_NAME)
    else:
        # Ensure index directory exists.
        os.makedirs(settings.INDEX_DIR, exist_ok=True)
        ix = create_in(dirname=settings.INDEX_DIR, schema=get_schema(), indexname=settings.INDEX_NAME)

    return ix


def delete_existing(ix, writer, uid):
    """
    Delete an existing post from the index.
    """

    searcher = ix.searcher()
    parser = MultifieldParser(fieldnames=['uid'], schema=ix.schema).parse(uid)
    indexed = searcher.search(parser)

    if not indexed.is_empty():
        # Delete the post from the index
        writer.delete_by_term('uid', uid, searcher=searcher)


def index_posts(posts, reindex=False):
    """
    Create or update a search index of posts.
    """
    # Indexes that already exist will to be updated
    # instead of starting from scratch.
    updating_index = index_exists()

    ix = init_index()

    # The writer is asynchronous by default
    writer = AsyncWriter(ix)
    elapsed, progress = timer_func()
    stream = islice(zip(count(1), posts), None)

    # Loop through posts and add to index
    for i, post in stream:
        progress(i, msg="posts indexed")
        # Delete an existing post before reindexing it.
        if updating_index:
            delete_existing(ix=ix, writer=writer, uid=post.uid)
        # Index post
        add_index(post=post, writer=writer)

    elapsed, progress = timer_func()

    # Commit additions to index

    if reindex:
        # Re-index posts when committing.
        writer.commit(mergetype=writing.CLEAR)
    else:
        writer.commit()

    elapsed(f"Created/updated index for {len(posts)} posts.")

    # Update indexed field on posts.
    posts.update(indexed=True)


def query(q='', fields=['content'], **kwargs):
    """
    Query the indexed, looking for a match in the specified fields.
    """
    ix = init_index()
    searcher = ix.searcher()

    parser = MultifieldParser(fieldnames=fields, schema=ix.schema).parse(q)
    results = searcher.search(parser, limit=settings.SEARCH_LIMIT, **kwargs)
    # Allow larger fragments
    results.fragmenter.maxchars = 300
    # Show more context before and after
    results.fragmenter.surround = 50

    #searcher.close()

    return results


