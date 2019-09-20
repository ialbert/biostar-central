import logging
from biostar.forum.models import Post
from biostar.forum import search
from django.conf import settings

logger = logging.getLogger("biostar")


def update_index(*args):
    """
    Index 1000 posts every 3 minutes
    """
    # Get un-indexed posts
    posts = Post.objects.filter(indexed=False)[:settings.BATCH_INDEXING_SIZE]

    # Nothing to be done.
    if not posts:
        logger.info("No new posts found")
        return

    logger.info(f"Indexing {len(posts)} posts.")

    # Update indexed field on posts.
    Post.objects.filter(id__in=posts.values('id')).update(indexed=True)

    try:
        search.index_posts(posts=posts)
        logger.info(f"Updated search index with {len(posts)} posts.")
    except Exception as exc:
        logger.info(f'Error updating index: {exc}')
        Post.objects.filter(id__in=posts.values('id')).update(indexed=False)

    return


if __name__ == '__main__':
    update_index()
