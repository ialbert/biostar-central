import logging
import time
from datetime import timedelta
from biostar.forum.models import Post
from biostar.accounts.models import Profile, User
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

    def wrap(*args, **kwargs):
        last = time.time()
        res = func(*args, **kwargs)
        diff = round(time.time() - last, 1)
        print(f'...function={func.__name__}; time={diff}secs')
        return res

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


def column_list(cursor):
    """
    Return list of column names currently found in the cursor.
    """
    colnames = [col[0] for col in cursor.description]
    return colnames


@timer
def select_related_users(posts, cursor):
    """
    Select all author and last edit users for each post
    """
    cols = {k: i for i, k in enumerate(column_list(cursor))}

    # Get author and last edit user ids into a set.
    user_ids = [(p[cols['author_id']], p[cols['lastedit_user_id']]) for p in posts if p]
    user_ids = {x for y in user_ids for x in y}
    user_ids = tuple(user_ids)

    # Joins the user profile as well.
    cursor.execute(f"""
                    SELECT *
                    FROM users_user
                    RIGHT JOIN users_profile ON (users_user.id = users_profile.user_id)
                    WHERE users_user.id IN {user_ids}
                    """)

    users = cursor.fetchall()
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


@timer
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
    post_cols = column_list(cursor=cursor)

    # Preform a select_related query for the users in each post.
    users = select_related_users(posts=threads, cursor=cursor)
    # Get the column name for the users table.
    user_cols = column_list(cursor=cursor)

    # Collapse results into a single dict
    context = dict(threads=dict(column=post_cols, rows=threads),
                   users=dict(column=user_cols, rows=users))

    logger.info(f"Start date={sd}")
    logger.info(f"End date={ed}")
    logger.info(f"Number of posts= {len(threads)}")
    logger.info(f"Number of users= {len(users)}")

    return context


@timer
def update_posts(threads, users):
    """
    Update local database with posts in 'threads'.
    Creates the posts if it does not exist.
    """
    rows = threads['rows']
    column = threads['column']
    elapsed, progress = util.timer_func()

    for row in rows:
        row = {col: val for col, val in zip(column, row)}
        exists = Post.objects.filter(uid=row['id']).exists()
        parent = Post.objects.filter(uid=row['parent_id']).first()
        root = Post.objects.filter(uid=row['root_id']).first()
        # lastedit_user = User.objects.filter(profile__uid=row['lastedit_user_id']).first()
        # author = User.objects.filter(profile__uid=row['author_id']).first()
        #
        # print(author, lastedit_user)
        #
        if not exists:

            post = Post.objects.create(uid=row['id'], creation_date=row['creation_date'],
                                       root=root, parent=parent,
                                       lastedit_user=users[row['lastedit_user_id']],
                                       author=users[row['author_id']],

                                       view_count=row['view_count'],
                                       vote_count=row['vote_count'],
                                       book_count=row['book_count'],
                                       accept_count=int(row['has_accepted']),
                                       thread_votecount=row['thread_score'],
                                       html=row['html'],
                                       content=row['content'],
                                       title=row['title'],
                                       status=row['status'],
                                       type=row['type'],
                                       tag_val=row['tag_val'],
                                       lastedit_date=row['lastedit_date'],
                                       rank=row['rank'])

    logger.info(f"Updated {len(rows)} posts.")
    return


@timer
def update_users(users):
    # Get the column names
    column = users.get('column')
    rows = users.get('rows', [])
    added = dict()
    # Check if this user already exists.
    for row in rows:
        # Map column names to row.
        row = {col: val for col, val in zip(column, row)}
        # Skip if user exists.
        if User.objects.filter(email=row['email']).exists():
            continue

        # Create the user
        username = f"{row['name'].replace(' ', '-')}-{row['user_id']}"
        user = User.objects.create(username=username, email=row['email'],
                                   password=row['password'], is_active=row['is_active'],
                                   is_staff=row['is_staff'], is_superuser=row['is_admin'])
        text = util.strip_tags(row['info'])
        # Update the Profile
        Profile.objects.filter(user=user).update(digest_prefs=row['digest_prefs'],
                                                 watched_tags=row['watched_tags'],
                                                 twitter=row['twitter_id'],
                                                 uid=row['user_id'], name=row['name'],
                                                 message_prefs=row['message_prefs'],
                                                 role=row['type'], last_login=row['last_login'],
                                                 html=row['info'], date_joined=row['date_joined'],
                                                 location=row['location'], website=row['website'],
                                                 scholar=row['scholar'], text=text,
                                                 score=row['score'], my_tags=row['my_tags'],
                                                 new_messages=row['new_messages'])
        added[row['user_id']] = user

    logger.info(f"Updated {len(rows)} users.")
    return added


@psycopg_required
@timer
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

    added = update_users(users=users)

    update_posts(threads=threads, users=added)

    return
