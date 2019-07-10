import logging
import copy
from django.conf import settings
from whoosh.qparser import MultifieldParser, OrGroup
from whoosh.analysis import SpaceSeparatedTokenizer, StopFilter, STOP_WORDS
from whoosh.index import create_in, open_dir
from whoosh.fields import ID, NGRAM, TEXT, KEYWORD, Schema, BOOLEAN, NUMERIC, NGRAMWORDS

from .models import Post
logger = logging.getLogger('engine')


# Stop words.
STOP = ['there', 'where', 'who'] + [w for w in STOP_WORDS]
STOP = set(STOP)


def add_index(post, writer):
    title = '' if not post.is_toplevel else post.title
    writer.add_document(title=title, url=post.get_absolute_url(),
                        type=post.get_type_display(),
                        lastedit_date=post.lastedit_date.timestamp(),
                        content=post.content, tags=post.tag_val,
                        is_toplevel=post.is_toplevel,
                        rank=post.rank, author=post.author.profile.name, uid=post.uid,
                        author_uid=post.author.profile.uid,
                        author_url=post.author.profile.get_absolute_url())


def open_index(index_dir=settings.INDEX_DIR, index_name=settings.INDEX_NAME):
    try:
        ix = open_dir(dirname=index_dir, indexname=index_name)
    except Exception as exc:
        logger.error(f"Error opening search index: {exc}")
        ix = None
    return ix


def create_index(posts, index_dir=settings.INDEX_DIR, index_name=settings.INDEX_NAME):
    """
    Create or update a search index of posts.
    """
    # Create the schema
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

    ix = create_in(index_dir, schema, indexname=index_name)
    # Exclude deleted posts from being indexed.
    posts = posts.exclude(status=Post.DELETED)
    writer = ix.writer()
    # Add the post info to index and commit.
    for post in posts:
        if not post.root:
            continue
        add_index(post=post, writer=writer)

    writer.commit()

    logger.info(f"Created index with: dir={index_dir}, name={index_name}")


def update_index(posts, index_dir=settings.INDEX_DIR, index_name=settings.INDEX_NAME):
    """Update indexed posts with the posts info"""
    # Open the index file we are about to update
    ix = open_index(index_dir=index_dir, index_name=index_name)

    # Prepare the searcher and writer objects
    searcher = ix.searcher()
    writer = ix.writer()

    for post in posts:
        # Check to see if this post already has an index
        parser = MultifieldParser(fieldnames=['uid'], schema=ix.schema).parse(post.uid)
        indexed = searcher.search(parser)
        # Delete already existing index
        if not indexed.is_empty():
            docnum = list(indexed.items())[0][0]
            writer.delete_document(docnum=docnum)

        # Reindex post
        add_index(post=post, writer=writer)

    writer.commit()
    searcher.close()
    logger.info(f"Updated index: updated posts={len(posts)} dir={index_dir}, name={index_name}")


def query(q='', fields=['content'], index_dir=settings.INDEX_DIR, index_name=settings.INDEX_NAME,
          **kwargs):

    ix = open_index(index_dir=index_dir, index_name=index_name)

    # def get_results():
    #     # TODO: trying to close the search stream and return results without raising ReaderClosed Exception
    #     #  deep-copying results did not work.
    #     print("foo")
    #     searcher = ix.searcher()
    #     def res():
    #         nonlocal query
    #         query = str(query)
    #         #query = query.encode('latin1')
    #         #print(query, type(query), query.encode("latin1"))
    #         q = MultifieldParser(fieldnames=fields, schema=ix.schema).parse(query)
    #         results = searcher.search(q, limit=settings.SIMILAR_FEED_COUNT, **kwargs)
    #         return results
    #     searcher.close()
    #     return res
    searcher = ix.searcher()

    parser = MultifieldParser(fieldnames=fields, schema=ix.schema).parse(q)
    results = searcher.search(parser, limit=settings.SIMILAR_FEED_COUNT, **kwargs)
    # Allow larger fragments
    results.fragmenter.maxchars = 300
    # Show more context before and after
    results.fragmenter.surround = 50

    #searcher.close()

    return results


