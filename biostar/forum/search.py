import logging
import time
import os
from itertools import count, islice

from django.conf import settings
from django.template import loader

from whoosh import writing
from whoosh.analysis import StemmingAnalyzer
from whoosh.sorting import FieldFacet, ScoreFacet
from whoosh.writing import AsyncWriter
from whoosh.qparser import MultifieldParser, OrGroup
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

    def progress(index, step=500, total=0, msg=""):
        nonlocal last
        if index % step == 0:
            percent = int((index / total) * 100) if total > index else index
            elapsed(f"... {percent}% ({index} out of {total}). {step} {msg}")

    return elapsed, progress


def index_exists():
    return exists_in(dirname=settings.INDEX_DIR, indexname=settings.INDEX_NAME)


def add_index(post, writer):
    writer.update_document(title=post.title, url=post.get_absolute_url(),
                           type_display=post.get_type_display(),
                           content_length=len(post.content),
                           type=post.type,
                           lastedit_date=post.lastedit_date.timestamp(),
                           content=post.content, tags=post.tag_val,
                           is_toplevel=post.is_toplevel,
                           rank=post.rank, uid=post.uid,
                           author_handle=post.author.username,
                           author=post.author.profile.name,
                           author_score=post.author.profile.score,
                           thread_votecount=post.thread_votecount,
                           vote_count=post.vote_count,
                           author_uid=post.author.profile.uid,
                           author_url=post.author.profile.get_absolute_url())


def get_schema():
    analyzer = StemmingAnalyzer(stoplist=STOP)
    schema = Schema(title=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    url=ID(stored=True),
                    content_length=NUMERIC(stored=True, sortable=True),
                    thread_votecount=NUMERIC(stored=True, sortable=True),
                    vote_count=NUMERIC(stored=True, sortable=True),
                    content=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    tags=KEYWORD(stored=True, commas=True),
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
                    type_display=TEXT(stored=True))
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


def index_posts(posts, overwrite=False):
    """
    Create or update a search index of posts.
    """

    ix = init_index()
    # The writer is asynchronous by default
    writer = AsyncWriter(ix)

    elapsed, progress = timer_func()
    total = posts.count()
    stream = islice(zip(count(1), posts), None)

    # Loop through posts and add to index
    for step, post in stream:
        progress(step, total=total, msg="posts indexed")
        add_index(post=post, writer=writer)

    # Commit to index
    if overwrite:
        logger.info("Overwriting the old index")
        writer.commit(mergetype=writing.CLEAR)
    else:
        writer.commit()

    elapsed(f"Indexed posts={total} dir={settings.INDEX_DIR} name={settings.INDEX_NAME}.")


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
    thread = FieldFacet('thread_votecount')
    content_length = FieldFacet("content_length", reverse=True)
    rank = FieldFacet("rank", reverse=True)
    default = ScoreFacet()

    # Splits the query into words and applies
    # and OR filter, eg. 'foo bar' == 'foo OR bar'
    orgroup = OrGroup

    # Sort by: toplevel, match score, author reputation, post rank.
    #sort_by = [post_type,  profile_score, rank, default]

    #sort_by = [post_type]

    #sort_by = [profile_score]

    #sort_by = [rank]

    #sort_by = [thread]

    sort_by = [default, content_length]

    #sort_by = [content_length]

    parser = MultifieldParser(fieldnames=fields, schema=ix.schema, group=orgroup).parse(q)
    results = searcher.search(parser, sortedby=sort_by, limit=settings.SEARCH_LIMIT, terms=True, **kwargs)
    # Allow larger fragments
    results.fragmenter.maxchars = 100
    #results.fragmenter.charlimit = None
    # Show more context before and after
    results.fragmenter.surround = 100

    return results
