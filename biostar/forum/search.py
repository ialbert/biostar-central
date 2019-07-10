import logging
import time
import os
from itertools import count, islice

from django.conf import settings
from whoosh.qparser import MultifieldParser, OrGroup
from whoosh.analysis import SpaceSeparatedTokenizer, StopFilter, STOP_WORDS
from whoosh.index import create_in, open_dir, exists_in
from whoosh.fields import ID, NGRAM, TEXT, KEYWORD, Schema, BOOLEAN, NUMERIC, NGRAMWORDS

from .models import Post
logger = logging.getLogger('engine')


# Stop words.
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


def init_index(mode='w'):

    os.makedirs(settings.INDEX_DIR, exist_ok=True)

    # Initialize a fresh index.
    if mode == 'w':
        ix = create_in(dirname=settings.INDEX_DIR, schema=get_schema(), indexname=settings.INDEX_NAME)

    index_exists = exists_in(dirname=settings.INDEX_DIR, indexname=settings.INDEX_DIR)

    if index_exists:
        # Open an existing index to update.
        ix = open_dir(dirname=settings.INDEX_DIR, indexname=settings.INDEX_NAME)
    else:
        ix = None
        pass

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
        docnum = list(indexed.items())[0][0]
        writer.delete_document(docnum=docnum)

    return


def index_posts(posts, ix=None):
    """
    Create or update a search index of posts.
    """

    if ix is None:
        ix = init_index()

    # Exclude deleted posts from being indexed.
    posts = posts.exclude(status=Post.DELETED)

    writer = ix.writer()
    elapsed, progress = timer_func()
    stream = islice(zip(count(1), posts), None)

    for i, post in stream:
        progress(i, msg="posts indexed")
        # Skip posts without a root,
        # happens when only transferring parts of the old biostar database.
        if not post.root:
            continue

        # Delete an existing post before reindexing it.
        if updating:
            delete_existing(ix=ix, writer=writer, uid=post.uid)
        # Index post
        add_index(post=post, writer=writer)

    elapsed, progress = timer_func()
    # Commit changes to the index.
    writer.commit()
    elapsed(f"Created/updated index for {len(posts)}")


def query(ix, q='', fields=['content'], **kwargs):
    """
    Query the indexed, looking for a match in the specified fields.
    """
    searcher = ix.searcher()

    parser = MultifieldParser(fieldnames=fields, schema=ix.schema).parse(q)
    results = searcher.search(parser, limit=settings.SIMILAR_FEED_COUNT, **kwargs)
    # Allow larger fragments
    results.fragmenter.maxchars = 300
    # Show more context before and after
    results.fragmenter.surround = 50

    #searcher.close()

    return results


