from django.conf import settings
from django.db.models import Q
from whoosh.writing import AsyncWriter
from whoosh.analysis import StemmingAnalyzer
from whoosh.fields import ID, TEXT, KEYWORD, Schema, NUMERIC
from biostar.forum.models import Post
from biostar.forum import search


def add_to_spam(post, writer):
    writer.update_document(title=post.title,
                           content_length=len(post.content),
                           content=post.content,
                           tags=post.tag_val,
                           uid=post.uid)


def spam_schema():
    analyzer = StemmingAnalyzer()
    schema = Schema(title=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    uid=ID(stored=True),
                    content_length=NUMERIC(stored=True, sortable=True),
                    content=TEXT(stored=True, analyzer=analyzer, sortable=True),
                    tags=KEYWORD(stored=True, commas=True))
    return schema


def build_spam_index(overwrite=False):

    # Get all un-indexed spam posts.
    spam = Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED), indexed=False)
    schema = spam_schema()

    # Initialize the spam index
    ix = search.init_index(dirname=settings.SPAM_INDEX_DIR, indexname=settings.SPAM_INDEX_NAME, schema=schema)

    # Index spam posts.
    search.index_posts(posts=spam, ix=ix, overwrite=overwrite, add_func=add_to_spam)

    return ix


def spam_score(post):
    """
    Calculate the spam score for a post.
    """

    # Get the average score for all spam hits found for this post.

    return


def search_spam(post):
    """
    Search spam index for posts similar to this one.
    Returns
    """

    # Build spam index or return already existing spam index.
    ix = build_spam_index()

    # Remove post from index once spam score is calculated.
    writer = AsyncWriter(ix)

    # Write this post to the index to perform most_like_this()
    add_to_spam(post=post, writer=writer)
    writer.commit()

    # Search for this post in the spam index
    fields = ['uid']
    results = search.preform_whoosh_search(ix=ix, query=post.uid, fields=fields)

    # Preform more_like_this on this posts content
    similar = results[0].more_like_this('content', top=500)

    # Remove this post from the spam index.
    writer = AsyncWriter(ix)
    writer.delete_by_term('uid', text=post.uid)
    writer.commit()

    # Get the results into a list
    similar = list(map(search.normalize_result, similar))
    results.searcher.close()
    return similar


def quarantine(post):
    """
    User + post pass a series of tests and rules to determine
    if they might be spam or not
    """

    # User's with high enough score automatically given green light.

    print(search_spam(post=post))

    # Mark this post as "maybe" being spam
    #Post.objects.filter(id=post.id).update(spam=Post.MAYBE_SPAM)













