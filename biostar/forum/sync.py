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


def timer(func):
    """
    Time a given function.
    """

    def wrap():


        return

    return


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


def get_column_names(cursor):
    """
    Return a dictionary of column names currently in the cursor.
    column name : column index
    """
    colnames = [col[0] for col in cursor.description]

    cols = {k: i for i, k in enumerate(colnames)}

    return cols


def select_related_users(posts, cursor):
    """
    Select all author and last edit users for each post
    """
    cols = get_column_names(cursor)

    # Get author and last edit user ids into a set.
    user_ids = [(p[cols['author_id']], p[cols['lastedit_user_id']]) for p in posts if p]
    user_ids = {x for y in user_ids for x in y}
    user_ids = tuple(user_ids)

    # Joins the user profile as well.
    cursor.execute(f"""
                    SELECT *
                    FROM users_user
                    INNER JOIN users_profile ON (users_user.id = users_profile.user_id)
                    WHERE users_user.id IN {user_ids}
                    """)

    users = cursor.fetchall()
    #colnames = [col[0] for col in cursor.description]
    #print([(idx, cols) for idx, cols in enumerate(colnames)])
    #1/0
    #cols = {k: i for i, k in enumerate(colnames)}
    #print(get_column_names(cursor), "KPPPPPPP")
    #1/0
    return users


def posts_query_str(start, end, batch=None):

    # Prefetch the thread associated with each post
    # Order by id so the roots and parents always come before children.
    q = f"""SELECT DISTINCT ON (thread.id) thread.id, thread.root_id, thread.parent_id, 
                                           thread.author_id, thread.lastedit_user_id, 
                                           thread.rank, thread.status, thread.type, 
                                           thread.creation_date, thread.lastedit_date,
                                           thread.title, thread.tag_val, 
                                           thread.content,thread.html, 
                                           thread.view_count, thread.vote_count, 
                                           thread.book_count, thread.has_accepted, 
                                           thread.thread_score
            FROM posts_post post
            INNER JOIN posts_post thread ON (post.root_id = thread.root_id)
            WHERE post.lastedit_date between '{start}'::timestamp and '{end}'::timestamp 
            OR post.creation_date between '{start}'::timestamp and '{end}'::timestamp  
            ORDER BY thread.id ASC
         """
    if batch:
        q += f"LIMIT '{batch}'"

    return q


def retrieve(cursor, start, days, batch=None):
    """
    Retrieve relevant data between a given date range.
    """
    # The end date is calculated
    end = start + timedelta(days=days)

    start, end = sorted([start, end])

    sd, ed = start.date(), end.date()

    query = posts_query_str(start=start, end=end, batch=batch)

    # Execute actual query and hit the database
    cursor.execute(query)
    threads = cursor.fetchall()

    # No posts exist for the given date range
    if not threads:
        logger.info(f"No posts found for start date={sd}, end date={ed}")
        return dict()

    # Retrieve column names from cursor to later map to rows
    post_cols = get_column_names(cursor=cursor)

    # Preform a select_related query for the users in each post.
    users = select_related_users(posts=threads, cursor=cursor)
    # Get the column name for the users table.
    user_cols = get_column_names(cursor=cursor)

    # Collapse results into a single dict
    context = dict(threads=dict(columns=post_cols, rows=threads),
                   users=dict(columns=user_cols, rows=users))

    logger.info(f"Start date={sd}")
    logger.info(f"End date={ed}")
    logger.info(f"Number of posts= {len(threads)}")
    logger.info(f"Number of users= {len(users)}")

    return context


def update_posts(threads):
    len(threads)

    # Get all of the users
    # print(threads[1])
    logger.info(f"Loading {len(threads)} posts.")

    # print(cur.fetchone())
    1 / 0

    # Get the user for this post, get the

    # Get all posts for a given root posts and load all posts for that root.

    # Load all root posts for these posts first.

    # Load non root posts afterwards.

    return


def update_users(users):

    # Get the column for the table
    column = users.get('columns')
    rows = users.get('rows', [])

    # Check if this user already exists.
    for row in rows:

        #user_dict = {name: val for name, val in zip(column, row)}
        print(row, sorted(column.items(), key=lambda x: x[1] ))

        print("FOOOOO" *5)
        1/0

    1/0

    return


@psycopg_required
def begin_sync(start=None, days=1, options=dict()):
    # Create initial connection to database
    conn = psycopg2.connect(dbname=options['dbname'],
                            host=options['host'],
                            user=options['user'],
                            password=options['password'],
                            port=options['port'],
                            sslmode='require')

    # Get the start date
    start = start or get_local_start_date(day_range=days)

    # Get the cursor.
    cur = conn.cursor()

    # Get all relevant data within a given timespan
    # data is a dict with threads, users, etc.
    data = retrieve(cursor=cur, start=start, days=days)

    users = data.get('users', {})
    threads = data.get('threads', {})

    for p in threads['rows']:
        print('-'*5)
        print("thread root, parent, post, post parent")
        print(p[0], p[1], p[2], p[3], p[4], )
        print('-'*5)

    #update_users(users=users)

    update_posts(threads=threads)

    return
