import os
import logging
from django.conf import settings
from whoosh.qparser import MultifieldParser, QueryParser
from whoosh.index import create_in, open_dir
from whoosh.fields import ID, NGRAM, TEXT, KEYWORD, Schema

logger = logging.getLogger('engine')


def create_index(posts, index_dir=settings.INDEX_DIR):

    schema = Schema(title=TEXT(stored=True), url=ID(stored=True), content=TEXT(stored=True),
                    tags=KEYWORD(stored=True))

    ix = create_in(index_dir, schema, indexname=settings.INDEX_NAME)
    writer = ix.writer()
    for post in posts:
        writer.add_document(title=f"{post.title}", url=post.get_absolute_url(),
                            content=post.content, tags=post.tag_val)
    writer.commit()

    logger.info(f"Created search index in {index_dir}")


def search_index(query='', index_dir=settings.INDEX_DIR):

    ix = open_dir(dirname=index_dir, indexname=settings.INDEX_NAME)

    with ix.searcher() as searcher:

        query = MultifieldParser(['content', 'title'], ix.schema).parse("This is a blog post")
        results = searcher.search(query)
        print([x for x in results])
        #if len(results):
        #    print(results[0].more_like_this(fieldname="content"))
        #    1/0

        return results




