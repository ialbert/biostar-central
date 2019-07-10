import logging
import copy
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


def index_posts(posts, clear=False, index_dir=settings.INDEX_DIR, index_name=settings.INDEX_NAME):
    """
    Create or update a search index of posts.
    """

    index_exists = exists_in(dirname=index_dir, indexname=index_name)
    if index_exists and not clear:
        # Open an existing index
        ix, updating = open_index(index_dir=index_dir, index_name=index_name), True
    else:
        # Create a brand new index
        ix, updating = create_in(dirname=index_dir, schema=get_schema(), indexname=index_name), False

    # Exclude deleted posts from being indexed.
    posts = posts.exclude(status=Post.DELETED)
    writer = ix.writer()

    for post in posts:
        # Skip posts without a root,
        # happens when only transferring parts of the biostar database.
        if not post.root:
            continue
        # Delete an existing post before reindexing it.
        if updating:
            delete_existing(ix=ix, writer=writer, uid=post.uid)
        # Index post
        add_index(post=post, writer=writer)

    # Commit changes to the index.
    writer.commit()
    logger.info(f"Created index with: dir={index_dir}, name={index_name}")


def query(q='', fields=['content'], index_dir=settings.INDEX_DIR, index_name=settings.INDEX_NAME,
          **kwargs):

    ix = open_index(index_dir=index_dir, index_name=index_name)

    searcher = ix.searcher()

    parser = MultifieldParser(fieldnames=fields, schema=ix.schema).parse(q)
    results = searcher.search(parser, limit=settings.SIMILAR_FEED_COUNT, **kwargs)
    # Allow larger fragments
    results.fragmenter.maxchars = 300
    # Show more context before and after
    results.fragmenter.surround = 50

    #searcher.close()

    return results


