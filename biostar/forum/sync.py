import logging
import time
import os
from datetime import timedelta, datetime
from functools import partial
from django.conf import settings
from django.db.models import Q
from biostar.forum.models import Post, Vote
from biostar.accounts.models import Profile, User, Logger
from biostar.forum import util
from django.core.cache import cache

try:
    import psycopg2

    PSYCOPG_INSTALLED = True
except ImportError as exc:
    PSYCOPG_INSTALLED = False

logger = logging.getLogger('engine')

# Most recent start date stored in a file for the next iteration
# TODO: being refactored out
START = os.path.join(settings.BASE_DIR, "export", 'start_sync.txt')


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

    def time_func(*args, **kwargs):
        last = time.time()
        res = func(*args, **kwargs)
        diff = round(time.time() - last, 1)
        logger.info(f"{func.__name__}() time = {diff}secs")
        return res

    return time_func


def column_list(cursor):
    """
    Return list of column names currently found in the cursor.
    """
    colnames = [col[0] for col in cursor.description]
    return colnames


@timer
def select_related_votes(posts, cursor):
    """
    Select votes that been made recently.
    """
    cols = {k: i for i, k in enumerate(column_list(cursor))}

    # Get author and last edit user ids into a set.
    post_ids = {p[cols['id']] for p in posts if p}
    post_ids = tuple(post_ids)

    cursor.execute(f"""
                    SELECT *
                    FROM posts_vote
                    INNER JOIN users_profile ON (posts_vote.author_id = users_profile.user_id)
                    WHERE posts_vote.post_id IN {post_ids}
                    """)

    votes = cursor.fetchall()
    return votes


def select_related_users(posts, cursor):
    """
    Select all author and last edit users found in list of posts.
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
        q += f" LIMIT '{batch}'"

    return q


def get_start():
    """
    Return the start date stored in cache.
    """
    if os.path.isfile(START):
        start = open(START, 'r').readline().strip()
        start = datetime.fromisoformat(start)
    else:
        start = util.now()

    return start


def set_start(start, end, days):

    # Start from the minimum of two time points when going backwards.
    if days < 0:
        store = min([start, end])
    else:
        store = max([start, end])

    open(START, 'w').write(str(store))


@timer
def retrieve(cursor, start, days, batch=None):
    """
    Retrieve relevant data between a given date range.
    """
    # The calculated end date
    end = start + timedelta(days=days)

    start, end = sorted([start, end])

    # Set the start day cache
    set_start(start, end, days)

    query = posts_query_str(start=start, end=end, batch=batch)

    # Execute query and hit the database
    cursor.execute(query)
    posts = cursor.fetchall()

    # No posts exist for the given date range
    if not posts:
        logger.info(f"No posts found for start={start.date()} end={end.date()}")
        return dict()

    # Retrieve column names from cursor to later map to rows
    post_cols = column_list(cursor=cursor)

    # Preform a select_related query for the users in each post.
    users = select_related_users(posts=posts, cursor=cursor)
    # Get the column name for the users table.
    user_cols = column_list(cursor=cursor)

    # Collapse results into a single dict
    context = dict(threads=dict(column=post_cols, rows=posts),
                   users=dict(column=user_cols, rows=users))

    threads = [r for r in posts if r[0] == r[1]]
    logger.info(f"Start\t{start.date()}")
    logger.info(f"End\t{end.date()}")
    logger.info(f"Number of total posts \t{len(posts)}")
    logger.info(f"Number of threads \t{len(threads)}")
    logger.info(f"Number of users \t{len(users)}")

    return context


def split_rows(rows):
    """
    Split incoming row into what is to be updated/created.
    """

    posts = set(Post.objects.all().values_list('uid', flat=True))
    update = list(filter(lambda r: str(r[0]) in posts, rows))
    create = list(filter(lambda r: str(r[0]) not in posts, rows))

    logger.info(f"Number being created\t{len(create)}")
    logger.info(f"Number being updated\t{len(update)}")
    return create, update


def bulk_create(rows, column, users, relations=dict()):
    for row in rows:
        # Map column names to row.
        row = {col: val for col, val in zip(column, row)}
        post = Post(lastedit_user=users[row['lastedit_user_id']],
                    author=users[row['author_id']],
                    uid=row['id'],
                    is_toplevel=row['id'] == row['root_id'],
                    view_count=row['view_count'],
                    vote_count=row['vote_count'],
                    book_count=row['book_count'],
                    accept_count=int(row['has_accepted']),
                    thread_votecount=row['thread_score'],
                    creation_date=row['creation_date'],
                    html=row['html'],
                    content=row['content'],
                    title=row['title'],
                    status=row['status'],
                    type=row['type'],
                    tag_val=row['tag_val'],
                    lastedit_date=row['lastedit_date'],
                    rank=row['rank'])
        relations[str(row['id'])] = str(row['root_id']), str(row['parent_id'])
        yield post


def set_relations(relations):

    posts = {post.uid: post for post in Post.objects.filter(uid__in=relations.keys())}
    for pid in relations:
        root_uid, parent_uid = relations[pid][0], relations[pid][1]
        post = posts[pid]
        root = posts.get(root_uid)
        parent = posts.get(parent_uid)
        if not (root and parent):
            print("MISSING", pid, root_uid, root)
            continue
        post.root = root
        post.parent = parent
        yield post


def set_counts(relations):

    descendants = lambda p: (Post.objects.filter(root=p).exclude(id=p.id)
                             if p.is_toplevel else Post.objects.filter(parent=p))

    posts = {post: (descendants(post).filter(type=Post.ANSWER).count(),
                    descendants(post).filter(type=Post.COMMENT).count())

             for post in Post.objects.filter(uid__in=relations.keys())}

    for post in posts:
        answer_count, comment_count = posts[post][0], posts[post][1]
        post.reply_count = answer_count + comment_count
        post.answer_count = answer_count
        post.comment_count = comment_count
        yield post


def bulk_create_votes(rows, pdict, udict):
    for row in rows:
        pass

def split_votes():

    return


@timer
def update_votes(rows, threads, users):

    posts_map = {post.uid: post for post in Post.objects.all()}
    users_map = {user.email: user for user in User.objects.all()}

    # Bulk create the votes.
    #Vote.objects.bulk_create(objs=bulk_create_votes(rows=rows, pdict, udict), batch_size=1000)

    return


@timer
def update_posts(threads, users, preform_updates=True):
    """
    Update local database with posts in 'threads'.
    Creates the posts if it does not exist.
    """
    rows = threads.get('rows', [])
    column = threads.get('column', [])
    relations = dict()

    # Split incoming rows into the
    to_create, to_update = split_rows(rows)

    # Bulk create new posts
    Post.objects.bulk_create(objs=bulk_create(rows=to_create, column=column, users=users,
                                              relations=relations),
                             batch_size=500)

    # # Bulk update existing posts
    # if preform_updates:
    #     #Post.objects.bulk_update(objs=updator, batch_size=500)
    #     pass
    # else:
    #     logger.info("Skipped preforming updates.")

    # Update the root, parent relationships in the queried posts
    Post.objects.bulk_update(objs=set_relations(relations=relations),
                             fields=["root", "parent"],
                             batch_size=1000)
    # Update counts on  posts
    Post.objects.bulk_update(objs=set_counts(relations=relations),
                             fields=["reply_count", "comment_count", "answer_count"],
                             batch_size=1000)

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
        user = User.objects.filter(email=row['email']).first()
        if user:
            added[row['user_id']] = user
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
def sync_db(start=None, days=1, options=dict()):
    # Create initial connection to database
    conn = psycopg2.connect(dbname=options['dbname'],
                            host=options['host'],
                            user=options['user'],
                            password=options['password'],
                            port=options['port'],
                            sslmode='require')

    # Get the start date from input or cache.
    # If none are provided, now() is returned.
    start = start or get_start()

    update = options['update']
    # Get the cursor.
    cur = conn.cursor()

    # Get all relevant data within a given timespan
    # data is a dict with threads, users, etc.
    data = retrieve(cursor=cur, start=start, days=days)

    users = data.get('users', {})
    threads = data.get('threads', {})
    votes = data.get('votes', {})

    added = update_users(users=users)

    update_posts(threads=threads, users=added, preform_updates=update)

    # update_votes()

    return
