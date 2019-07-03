import os
import copy
import logging
from django.conf import settings
from whoosh.qparser import MultifieldParser, QueryParser, OrGroup
from whoosh.index import create_in, open_dir
from whoosh.fields import ID, NGRAM, TEXT, KEYWORD, Schema, BOOLEAN, NUMERIC

logger = logging.getLogger('engine')


def create_index(posts, index_dir=settings.INDEX_DIR, index_name=settings.INDEX_NAME):
    """
    Create search index of posts.
    Created inside of directory with the
    """
    # Create the schema
    schema = Schema(title=TEXT(stored=True), url=ID(stored=True), content=TEXT(stored=True),
                    tags=KEYWORD(stored=True), toplevel=BOOLEAN(stored=True),
                    rank=NUMERIC(stored=True, sortable=True), author=TEXT(stored=True),
                    author_url=ID(stored=True))

    ix = create_in(index_dir, schema, indexname=index_name)
    writer = ix.writer()
    # Add the post contents to index and commit.
    for post in posts:
        writer.add_document(title=f"{post.title}", url=post.get_absolute_url(),
                            content=post.content, tags=post.tag_val, toplevel=post.is_toplevel,
                            rank=post.rank, author=post.author.profile.name,
                            author_url=post.author.profile.get_absolute_url())
    writer.commit()

    logger.info(f"Created search index in {index_dir}")


def search_index(query='', fields=['content'], index_dir=settings.INDEX_DIR):

    ix = open_dir(dirname=index_dir, indexname=settings.INDEX_NAME)

    #with ix.searcher() as searcher:
    searcher = ix.searcher()
    # Group each word with an OR like such:
    # with a query 'foo bar' search 'foo OR bar'
    og = OrGroup
    query = MultifieldParser(fields, ix.schema, group=og).parse(query)
    results = searcher.search(query, sortedby='rank', reverse=True)
    #searcher.close()
    return results




