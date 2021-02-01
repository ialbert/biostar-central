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
from whoosh.writing import AsyncWriter
from whoosh.searching import Results

import html2markdown
import bleach
from whoosh.qparser import MultifieldParser, OrGroup
from whoosh.analysis import STOP_WORDS
from whoosh.index import create_in, open_dir, exists_in
from whoosh.fields import ID, TEXT, KEYWORD, Schema, BOOLEAN, NUMERIC, DATETIME

from biostar.forum.models import Post

logger = logging.getLogger('biostar')

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
        print(f"{msg} in {sec} seconds")

    def progress(index, step=500, total=0, msg=""):
        nonlocal last
        if index % step == 0:
            percent = int((index / total) * 100) if total >= index else index
            elapsed(f"... {percent}% ({index} out of {total}). {step} {msg}")

    return elapsed, progress


class SearchResult(object):
    def __init__(self, **kwargs):
        self.title = ''
        self.content = ''
        self.total = 0
        self.url = ''
        self.type_display = ''
        self.content_length = ''
        self.type = ''
        self.lastedit_date = ''
        self.is_toplevel = ''
        self.rank = ''
        self.pagenum = 0
        self.pagecount = 0
        self.uid = ''
        self.author_handle = ''
        self.author = ''
        self.author_email = ''
        self.author_score = ''
        self.thread_votecount = ''
        self.vote_count = ''
        self.author_uid = ''
        self.author_url = ''
        self.__dict__.update(kwargs)

    def __iter__(self):
        yield

    def is_last_page(self):
        return True

    def __len__(self):
        return self.total


def normalize_result(result):
    "Return a bunch object for result."

    # Result is a database object.
    bunched = SearchResult(title=result.get('title'), content=result.get('content'), url=result.get('url'),
                           type_display=result.get('type_display'), content_length=result.get('content_length'),
                           type=result.get('type'), lastedit_date=result.get('lastedit_date'),
                           is_spam=result.get("is_spam", False), is_toplevel=result.get('is_toplevel'),
                           rank=result.get('rank'), uid=result.get('uid'),
                           author_handle=result.get('author_handle'), author=result.get('author'),
                           author_score=result.get('author_score'), score=result.score,
                           thread_votecount=result.get('thread_votecount'), vote_count=result.get('vote_count'),
                           author_uid=result.get('author_uid'), author_url=result.get('author_url'))

    return bunched


def index_exists(dirname=settings.INDEX_DIR, indexname=settings.INDEX_NAME):
    return exists_in(dirname=dirname, indexname=indexname)


def add_index(post, writer):

    # Ensure the content is stripped of any html.
    content = bleach.clean(post.content, styles=[], attributes={}, tags=[], strip=True)
    writer.update_document(title=post.title, url=post.get_absolute_url(),
                           type_display=post.get_type_display(),
                           content_length=len(content),
                           type=post.type,
                           creation_date=post.creation_date,
                           lastedit_date=post.lastedit_date,
                           lastedit_user=post.lastedit_user.profile.name,
                           lastedit_user_email=post.author.email,
                           lastedit_user_score=post.author.profile.score,
                           lastedit_user_uid=post.author.profile.uid,
                           lastedit_user_url=post.lastedit_user.profile.get_absolute_url(),
                           content=content,
                           tags=post.tag_val,
                           is_toplevel=post.is_toplevel,
                           rank=post.rank, uid=post.uid,
                           vote_count=post.vote_count,
                           reply_count=post.reply_count,
                           view_count=post.view_count,
                           author_handle=post.author.username,
                           author=post.author.profile.name,
                           answer_count=post.root.answer_count,
                           root_has_accepted=post.root.has_accepted,
                           author_email=post.author.email,
                           author_score=post.author.profile.score,
                           thread_votecount=post.thread_votecount,
                           author_uid=post.author.profile.uid,
                           author_url=post.author.profile.get_absolute_url(),
                           author_is_moderator=post.author.profile.is_moderator,
                           author_is_suspended=post.author.profile.is_suspended,
                           lastedit_user_is_suspended=post.lastedit_user.profile.is_suspended,
                           lastedit_user_is_moderator=post.lastedit_user.profile.is_moderator)


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
                    author_is_moderator=BOOLEAN(stored=True),
                    lastedit_user_is_moderator=BOOLEAN(stored=True),
                    lastedit_user_is_suspended=BOOLEAN(stored=True),
                    author_is_suspended=BOOLEAN(stored=True),
                    lastedit_date=DATETIME(stored=True, sortable=True),
                    creation_date=DATETIME(stored=True, sortable=True),
                    rank=NUMERIC(stored=True, sortable=True),
                    author=TEXT(stored=True),
                    lastedit_user=TEXT(stored=True),
                    lastedit_user_email=TEXT(stored=True),
                    lastedit_user_score=NUMERIC(stored=True, sortable=True),
                    lastedit_user_uid=ID(stored=True),
                    lastedit_user_url=ID(stored=True),
                    author_score=NUMERIC(stored=True, sortable=True),
                    author_handle=TEXT(stored=True),
                    author_email=TEXT(stored=True),
                    author_uid=ID(stored=True),
                    author_url=ID(stored=True),
                    root_has_accepted=BOOLEAN(stored=True),
                    reply_count=NUMERIC(stored=True, sortable=True),
                    view_count=NUMERIC(stored=True, sortable=True),
                    answer_count=NUMERIC(stored=True, sortable=True),
                    uid=ID(stored=True),
                    type=NUMERIC(stored=True, sortable=True),
                    type_display=TEXT(stored=True))
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


def print_info(dirname=None, indexname=None,):
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
        writer.commit()

    elapsed(f"Indexed posts={total}")


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


def preform_whoosh_search(query, ix=None, fields=None, page=None, per_page=None, sortedby=[], reverse=True,
                          **kwargs):
    """
        Query the indexed, looking for a match in the specified fields.
        Results a tuple of results and an open searcher object.
        """

    per_page = per_page or settings.SEARCH_RESULTS_PER_PAGE
    fields = fields or ['tags', 'title', 'author', 'author_uid', 'content', 'author_handle']
    ix = ix or init_index()
    searcher = ix.searcher()

    # Splits the query into words and applies
    # and OR filter, eg. 'foo bar' == 'foo OR bar'
    orgroup = OrGroup

    parser = MultifieldParser(fieldnames=fields, schema=ix.schema, group=orgroup).parse(query)
    if page:
        # Return a pagenated version of the results.
        results = searcher.search_page(parser,
                                       pagenum=page, pagelen=per_page, sortedby=sortedby,
                                       reverse=reverse,
                                       terms=True)
        results.results.fragmenter.maxchars = 100
        # Show more context before and after
        results.results.fragmenter.surround = 100
    else:
        results = searcher.search(parser, limit=settings.SEARCH_LIMIT, sortedby=sortedby, reverse=reverse,
                                  terms=True)
        # Allow larger fragments
        results.fragmenter.maxchars = 100
        results.fragmenter.surround = 100

    #logger.info("Preformed index search")

    return results


def preform_search(query, fields=None, top=0, sortedby=[], more_like_this=False):

    top = top or settings.SIMILAR_FEED_COUNT
    length = len(query.replace(" ", ""))

    if length < settings.SEARCH_CHAR_MIN:
        return []
    fields = fields or ['tags', 'title', 'author', 'author_uid', 'author_handle']
    whoosh_results = preform_whoosh_search(query=query, sortedby=sortedby, fields=fields)

    if more_like_this and len(whoosh_results):
            results = whoosh_results[0].more_like_this("content", top=top)
            # Filter results for toplevel posts.
            results = list(filter(lambda p: p['is_toplevel'] is True, results))
    else:
        results = whoosh_results

    # Ensure returned results types stay consistent.
    final_results = list(map(normalize_result, results))
    if not len(final_results):
        return []

    # Ensure searcher object gets closed.
    whoosh_results.searcher.close()

    return final_results
