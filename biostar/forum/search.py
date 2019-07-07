import os
import copy
import whoosh.query as search_query
import logging
from django.conf import settings
from whoosh.qparser import MultifieldParser, QueryParser, OrGroup
from whoosh.analysis import SpaceSeparatedTokenizer
from whoosh.index import create_in, open_dir
from whoosh.fields import ID, NGRAM, TEXT, KEYWORD, Schema, BOOLEAN, NUMERIC, NGRAMWORDS

from .models import Post
logger = logging.getLogger('engine')


def create_index(posts, index_dir=settings.INDEX_DIR, index_name=settings.INDEX_NAME):
    """
    Create search index of posts.
    Created inside of directory with the
    """
    # Create the schema
    tokenizer = SpaceSeparatedTokenizer()
    schema = Schema(title=NGRAMWORDS(stored=True, tokenizer=tokenizer), url=ID(stored=True),
                    content=NGRAMWORDS(stored=True, at='start', tokenizer=tokenizer),
                    tags=KEYWORD(stored=True), is_toplevel=BOOLEAN(stored=True), author_uid=ID(stored=True),
                    rank=NUMERIC(stored=True, sortable=True), author=TEXT(stored=True),
                    author_url=ID(stored=True), uid=ID(stored=True))

    ix = create_in(index_dir, schema, indexname=index_name)
    # Exclude deleted posts from being indexed.
    posts = posts.exclude(status=Post.DELETED)
    writer = ix.writer()
    # Add the post info to index and commit.
    for post in posts:
        if not post.root:
            continue
        writer.add_document(title=f"{post.title}", url=post.get_absolute_url(),
                            content=post.content, tags=post.tag_val, is_toplevel=post.is_toplevel,
                            rank=post.rank, author=post.author.profile.name, uid=post.uid,
                            author_uid=post.author.profile.uid, author_url=post.author.profile.get_absolute_url())
    writer.commit()

    logger.info(f"Created index with: dir={index_dir}, name={index_name}")


def search_index(query='', fields=['content'], index_dir=settings.INDEX_DIR, index_name=settings.INDEX_NAME, **kwargs):

    ix = open_dir(dirname=index_dir, indexname=index_name)

    #with ix.searcher() as searcher:
    searcher = ix.searcher()
    # Group each word with an OR clause.
    # For example, the query 'foo bar' will be: 'foo OR bar'
    og = OrGroup
    query = MultifieldParser(fields, ix.schema, group=og).parse(query)
    results = searcher.search(query, limit=settings.SIMILAR_FEED_COUNT, **kwargs)
    #searcher.close()
    return results




