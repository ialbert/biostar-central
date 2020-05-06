import logging
from datetime import timedelta
from biostar.forum.models import Post
from biostar.forum import util

try:
    import psycopg2

    PSYCOPG_INSTALLED = True
except ImportError as exc:
    PSYCOPG_INSTALLED = False

logger = logging.getLogger('engine')


def psycopg_required(func):
    """
    Ensure psycopg2 is installed before calling function
    """

    def wrap(*args, **kwagrs):
        if PSYCOPG_INSTALLED:
            return func(*args, **kwagrs)
        logger.error(f"psycopg2 not installed.")
        return

    return wrap


def get_local_start_date(day_range):
    """
    Return creation date of the most recently synced post in local database.
    """

    if day_range < 0:
        # Get the oldest post to continue syncing when day_range is negative
        post = Post.objects.old().order_by('creation_date').first()
    else:
        # Get the newest post to continue syncing when day_range is positive
        post = Post.objects.old().order_by('-creation_date').first()

    start_date = post.creation_date if post else util.now()
    return start_date


def query_posts(cursor, start, end):

    start, end = sorted([start, end])
    logger.info(f"Start date={start.date()}")
    logger.info(f"End date={end.date()}")

    # Gets all root and parent posts.
    cursor.execute(f"""SELECT *
                       FROM posts_post A
                       INNER JOIN posts_post B ON (A.root_id = B.id)
                       WHERE A.creation_date between '{start}'::timestamp and '{end}'::timestamp
                       """)

    return cursor.fetchall()


def get_closest_dates(cursor, start, days):

    logger.info(f"Getting the closest start date to: {start.date()}")
    end = start + timedelta(days=days)

    cursor.execute(f"""SELECT creation_date
                    FROM posts_post
                    WHERE creation_date < '{start}'::timestamp 
                    ORDER BY creation_date DESC
                    LIMIT 1;
                    """)

    new_start = cursor.fetchall()
    if not new_start:
        return start, end

    start = new_start[0][0]
    end = start + timedelta(days=days)

    return start, end


def load_posts(posts):
    print(posts[0])
    logger.info(f"Loading {len(posts)} posts.")

    # print(cur.fetchone())
    1 / 0

    # Get the user for this post, get the

    # Get all posts for a given root posts and load all posts for that root.

    # Load all root posts for these posts first.

    # Load non root posts afterwards.

    return


def get_posts(connection, start, days, batch=None):
    """
    """

    # Get the start and end date
    start = start or get_local_start_date(day_range=days)
    end = start + timedelta(days=days)

    # Get the cursor.
    cur = connection.cursor()

    # Get all of the posts for given time range
    posts = query_posts(cursor=cur, start=start, end=end)

    # There are no posts for the given time range.
    if not posts:
        # Get the closest times possible and retry.
        start, end = get_closest_dates(cursor=cur, start=start, days=days)

        # Query using new start and end times
        posts = query_posts(cursor=cur, start=start, end=end)

    return posts


@psycopg_required
def begin_sync(start=None, days=1, options=dict()):

    # Create initial connection to database
    conn = psycopg2.connect(dbname=options['dbname'],
                            host=options['host'],
                            user=options['user'],
                            password=options['password'],
                            port=options['port'],
                            sslmode='require')

    # Preform database query to get posts.
    posts = get_posts(connection=conn, start=start, days=days)

    load_posts(posts=posts)

    return
