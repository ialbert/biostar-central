import logging
import time
import os
from datetime import timedelta, datetime
from functools import partial
from django.conf import settings
from django.db.models import Q
from biostar.forum.models import Post, Vote, Sync, Subscription
from biostar.accounts.models import Profile, User, Logger
from biostar.forum import util, auth
from django.core.cache import cache

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


def select_related(posts, cursor, fr='posts_vote', where='posts_vote.post_id'):
    cols = {k: i for i, k in enumerate(column_list(cursor))}

    # Get author and last edit user ids into a set.
    post_ids = {p[cols['id']] for p in posts if p}
    post_ids = tuple(post_ids)

    cursor.execute(f"""
                    SELECT *
                    FROM {fr}
                    WHERE {where} IN {post_ids}
                    """)

    objs = cursor.fetchall()
    return objs


def select_related_votes(posts, cursor):
    """
    Select votes that been made recently.
    """
    from_q = 'posts_vote'
    where_q = 'posts_vote.post_id'
    return select_related(posts=posts, cursor=cursor, fr=from_q, where=where_q)


def select_related_subs(posts, cursor):
    from_q = 'posts_subscription'
    where_q = 'posts_subscription.post_id'
    return select_related(posts=posts, cursor=cursor, fr=from_q, where=where_q)


def select_related_users(user_ids, cursor):
    """
    Select all author and last edit users found in list of posts.
    """
    # Get author and last edit user ids into a set.
    user_ids = {str(x) for y in user_ids for x in y}
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


def most_recent(cursor):

    # Fetch the most current post made in the remote database
    cursor.execute(f"""
                     SELECT *
                     FROM posts_post
                     ORDER BY posts_post.id DESC 
                     LIMIT 1
                     """)

    post = cursor.fetchone()
    cols = {k: i for i, k in enumerate(column_list(cursor))}

    # Get the creation date for this post.
    date = post[cols['creation_date']]

    return date


def get_start(days, cursor):
    """
    Return the start date stored in database
    """

    if days < 0:
        recent = Sync.objects.filter(pk=1).first()
        recent = recent.last_synced if recent else most_recent(cursor=cursor)
        recent = recent + timedelta(days=2)
    else:
        # Get the most recently created post on the remote server.
        recent = Post.objects.old().order_by('-creation_date').first()
        recent = recent.creation_date if recent else None
        recent = recent - timedelta(days=2)

    start = recent or util.now()

    return start


def store_synced_date(date):
    """
    Stored when syncing backwards.
    """
    recent = Sync.objects.filter(pk=1).first()

    # Go back 12 hours from this date to ensure nothing is missed.
    date = date + timedelta(hours=12)

    # Update the existing sync date
    if recent:
        Sync.objects.filter(pk=1).update(last_synced=date)
    else:
        Sync.objects.create(last_synced=date)


def split_rows(rows):
    """
    Split incoming row into what is to be updated/created.
    """

    posts = set(Post.objects.all().values_list('uid', flat=True))
    update = list(filter(lambda r: str(r[0]) in posts, rows))
    create = list(filter(lambda r: str(r[0]) not in posts, rows))
    threads = list(filter(lambda r: (str(r[0]) not in posts) and r[0] == r[1], rows))

    logger.info(f"Number being of posts created\t{len(create)}")
    logger.info(f"Number being of threads created\t{len(threads)}")
    logger.info(f"Number being of posts updated\t{len(update)}")
    return create, update


def bulk_create_posts(rows, column, users, relations=dict()):
    for row in rows:
        # Map column names to row.
        row = {col: val for col, val in zip(column, row)}
        post = Post(lastedit_user=users[row['lastedit_user_id']],
                    author=users[row['author_id']],
                    uid=row['id'],
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


def bulk_relations(relations):
    posts = {post.uid: post for post in Post.objects.filter(uid__in=relations.keys())}
    for pid in relations:
        root_uid, parent_uid = relations[pid][0], relations[pid][1]
        post = posts[pid]
        root = posts.get(root_uid)
        parent = posts.get(parent_uid)
        if not (root and parent):
            continue
        post.root = root
        post.parent = parent
        post.is_toplevel = root == post
        yield post


def bulk_counts(relations):
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


def bulk_create_votes(rows, column, pdict, udict):
    for row in rows:
        row = {col: val for col, val in zip(column, row)}
        post = pdict[str(row['post_id'])]
        author = udict[row['author_id']]
        vtype = row['type']
        # Skip incomplete post/author information.
        if not (post and author) or not post.root:
            continue

        vote = Vote(post=post, author=author, type=vtype, uid=row['id'], date=row['date'])

        yield vote


def bulk_create_subs(rows, column, pdict, udict):
    for row in rows:
        row = {col: val for col, val in zip(column, row)}
        user = udict.get(row['user_id'])
        post = pdict.get(str(row['post_id']))

        # Skip incomplete data.
        if not (user and post):
            continue

        sub = Subscription(uid=row['id'], type=row['type'], user=user, post=post, date=row['date'])

        yield sub


def clean_votes(rows, column):
    """
    Delete votes that already exist in preparation for bulk creating.
    """
    for row in rows:
        row = {col: val for col, val in zip(column, row)}
        user_uid = row['author_id']
        post_uid = row['post_id']
        vtype = row['type']
        # Check if votes exists with the post and user.
        vote = Vote.objects.filter(author__profile__uid=user_uid, post__uid=post_uid, type=vtype).first()
        # Delete the existing vote and
        if vote:
            vote.delete()


def clean_subs(rows, column):
    """
    Delete subscriptions that already exist in preparation for bulk creating.
    """
    for row in rows:
        row = {col: val for col, val in zip(column, row)}
        user_uid = row['user_id']
        post_uid = row['post_id']
        # Check if sub exists with this post and user
        sub = Subscription.objects.filter(user__profile__uid=user_uid, post__uid=post_uid).first()
        # Delete the existing sub and
        if sub:
            sub.delete()


@timer
def sync_votes(votes, users):
    rows = votes.get('rows', [])
    column = votes.get('column', [])

    cols = {k: i for i, k in enumerate(column)}

    post_ids = set([str(r[cols.get('post_id')]) for r in rows])

    posts = {p.uid: p for p in Post.objects.filter(uid__in=post_ids)}

    # Clean the votes by deleting existing ones.
    clean_votes(rows=rows, column=column)

    # Bulk create votes.
    generator = bulk_create_votes(rows=rows, column=column, pdict=posts, udict=users)

    Vote.objects.bulk_create(objs=generator, batch_size=500)


@timer
def sync_subs(subs, users):
    rows = subs.get('rows', [])
    column = subs.get('column', [])

    cols = {k: i for i, k in enumerate(column)}

    post_ids = set([str(r[cols.get('post_id')]) for r in rows])

    posts = {p.uid: p for p in Post.objects.filter(uid__in=post_ids)}

    # Clean the subs by deleting existing ones.
    clean_subs(rows=rows, column=column)

    # Bulk create subs
    generator = bulk_create_subs(rows=rows, pdict=posts, column=column, udict=users)

    Subscription.objects.bulk_create(objs=generator, batch_size=500)


def get_user_ids(rows, column, is_posts=False, is_subs=False):
    # Get all user id's found in a a row.
    cols = {k: i for i, k in enumerate(column)}

    if is_posts:
        user_ids = [(r[cols['author_id']], r[cols['lastedit_user_id']]) for r in rows if r]
    elif is_subs:
        user_ids = [(r[cols['user_id']],) for r in rows if r]
    else:
        user_ids = [(r[cols['author_id']],) for r in rows if r]

    return user_ids


@timer
def retrieve(cursor, start, days, batch=None):
    """
    Retrieve relevant data between a given date range.
    """
    # The calculated end date
    end = start + timedelta(days=days)
    start, end = sorted([start, end])
    query = posts_query_str(start=start, end=end, batch=batch)

    # Execute query and return posts.
    cursor.execute(query)
    posts = cursor.fetchall()

    # No posts exist for the given date range
    if not posts:
        logger.info(f"No posts found for start={start.date()} end={end.date()}")
        return dict()

    # Retrieve column names from cursor to later map to rows
    post_cols = column_list(cursor=cursor)

    # Get the user ids involved with this post.
    user_ids = get_user_ids(rows=posts, is_posts=True, column=post_cols)

    # Preform a select_related query for the users and votes.
    votes = select_related_votes(posts=posts, cursor=cursor)
    votes_cols = column_list(cursor=cursor)

    # Add user id's involved with votes after the 'cursor' has been altered.
    user_ids += get_user_ids(rows=votes, column=votes_cols)

    # Preform a select_related query for the users and votes.
    subs = select_related_subs(posts=posts, cursor=cursor)
    subs_cols = column_list(cursor=cursor)

    user_ids += get_user_ids(rows=subs, is_subs=True, column=subs_cols)

    # Select all related users found in votes and subs
    users = select_related_users(user_ids=user_ids, cursor=cursor)
    user_cols = column_list(cursor=cursor)

    # Store the last synced date into the database.
    store_synced_date(start)

    # Collapse results into a single dict
    context = dict(threads=dict(column=post_cols, rows=posts),
                   users=dict(column=user_cols, rows=users),
                   votes=dict(column=votes_cols, rows=votes),
                   subs=dict(column=subs_cols, rows=subs))

    threads = [r for r in posts if r[0] == r[1]]
    logger.info(f"Start\t{start.date()}")
    logger.info(f"End\t{end.date()}")
    logger.info(f"Number of total posts \t{len(posts)}")
    logger.info(f"Number of threads \t{len(threads)}")
    logger.info(f"Number of votes \t{len(votes)}")
    logger.info(f"Number of subscriptions \t{len(subs)}")
    logger.info(f"Number of users \t{len(users)}")

    return context


@timer
def sync_posts(threads, users, update=True):
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
    creator_gen = bulk_create_posts(rows=to_create, column=column, users=users, relations=relations)
    Post.objects.bulk_create(objs=creator_gen, batch_size=500)

    # Bulk update existing posts
    if update:
        # Post.objects.bulk_update(objs=updator, batch_size=500)
        pass
    else:
        logger.info("Skipped updates.")

    # Update the root, parent relationships in the queried posts
    Post.objects.bulk_update(objs=bulk_relations(relations=relations),
                             fields=["root", "parent", 'is_toplevel'],
                             batch_size=1000)
    # Update counts on  posts
    Post.objects.bulk_update(objs=bulk_counts(relations=relations),
                             fields=["reply_count", "comment_count", "answer_count"],
                             batch_size=1000)


@timer
def sync_users(users):
    """
    Sync users one at a time instead of bulk creating.
    """
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
    # Get the cursor.
    cur = conn.cursor()

    # Get the start date from input or cache.
    # If none are provided, now() is returned.
    start = start or get_start(days=days, cursor=cur)

    update = options['update']

    # Get all relevant data within a given timespan
    # data is a dict with threads, users, etc.
    data = retrieve(cursor=cur, start=start, days=days)

    users = data.get('users', {})
    threads = data.get('threads', {})
    votes = data.get('votes', {})
    subs = data.get('subs', {})

    # Sync users fist
    synced_users = sync_users(users=users)

    # Then sync posts, votes, and subscriptions, and awards.
    sync_posts(threads, users=synced_users, update=update)

    sync_votes(votes, users=synced_users)

    sync_subs(subs, users=synced_users)

    return


@timer
def db_report(cursor, synced):

    # Get a count of all of the posts
    synced = tuple(synced)
    nposts = f"SELECT COUNT(*) FROM posts_post WHERE posts_post.id NOT IN {synced}"

    cursor.execute(nposts)
    nposts = cursor.fetchone()[0]

    # Return the newest post date.
    newest = f""" SELECT posts_post.creation_date FROM posts_post 
                       WHERE posts_post.id NOT IN {synced} 
                       ORDER BY posts_post.id DESC LIMIT 1"""
    cursor.execute(newest)
    newest = cursor.fetchone()[0]

    return nposts, newest


@psycopg_required
@timer
def report(start=None, days=1, options=dict()):
    # Create initial connection to database
    conn = psycopg2.connect(dbname=options['dbname'],
                            host=options['host'],
                            user=options['user'],
                            password=options['password'],
                            port=options['port'],
                            sslmode='require')

    # Get the cursor.
    cur = conn.cursor()

    # Get the start date from input or cache.
    # If none are provided, now() is returned.
    start = start or get_start(days=days, cursor=cur)

    retrieve(cursor=cur, start=start, days=days)

    already_synced = Post.objects.old().values_list('uid', flat=True)

    missing, newest = db_report(cursor=cur, synced=already_synced)
    # Get the currently loaded data
    logger.info(f"Synced posts\t{already_synced.count()}")
    logger.info(f"Missing posts \t{missing}")
    logger.info(f"Newest post \t{newest}")

    return
