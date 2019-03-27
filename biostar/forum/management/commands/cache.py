import requests
import hjson
import time
import logging
import os
from calendar import Calendar
from itertools import chain
from urllib.parse import urljoin

from django.core.management.base import BaseCommand
from django.conf import settings
from biostar.forum import util
from biostar.forum.models import Vote, Post
from biostar.accounts.models import User, Profile

logger = logging.getLogger(settings.LOGGER_NAME)


def get_data(full_url):
    """
    Fetch data from remote url
    """
    while True:

        try:
            # 5 sec timeout then retry
            response = requests.get(full_url, timeout=5)
            logger.info(f"Hit remote site:{full_url}")
            break
        except Exception as exc:
            logger.error(f"{exc}...sleeping for 5 seconds then retrying.")
            time.sleep(2)

    data = {} if response.status_code == 404 else hjson.loads(response.text)
    return data, response


def cache_post(root_dir, post_list, prefix):

    api_url = "https://www.biostars.org/api/post/"
    fname = f"{prefix}_posts.hjson"
    fullpath = os.path.join(root_dir, fname)
    cached_posts = hjson.loads(open(fullpath, 'r')) if os.path.exists(fullpath) else []

    for post_id in post_list:
        # Build url and fetch data for post
        full_url = urljoin(api_url, f"{post_id}")
        data, response = get_data(full_url=full_url)

        # No data found for given post id
        if not data:
            logger.warning(f"Post Id {post_id} does not exist.")
            continue
        # Cache posts.
        cached_posts.append(hjson.dumps(data))

    # Write to file.
    with open(fullpath, "w") as cache:
        hjson.dump(cached_posts, cache)
        logger.info(f"Cached {len(cached_posts)} posts into {fullpath}")

    return


def cache_users(root_dir, user_list, prefix):

    api_url = "https://www.biostars.org/api/user/"
    # Append to existing file.
    fname = f"{prefix}_users.hjson"
    fullpath = os.path.join(root_dir, fname)
    cached_users = hjson.load(open(fullpath, 'r')) if os.path.exists(fullpath) else []

    for user_id in user_list:
        # Build url and fetch data for post
        full_url = urljoin(api_url, f"{user_id}")
        data, response = get_data(full_url=full_url)

        # No data found for given post id
        if not data:
            logger.warning(f"User Id {user_id} does not exist.")
            continue
        # Add to cache
        cached_users.append(hjson.dumps(data))

    # Write to file.
    with open(fullpath, "w") as cache:
        hjson.dump(cached_users, cache)
        logger.info(f"Cached {len(cached_users)} users into {fullpath}")

    return


def cache_votes(root_dir, votes_list, prefix):

    api_url = "https://www.biostars.org/api/vote/"
    fname = f"{prefix}_votes.hjson"
    fullpath = os.path.join(root_dir, fname)
    cached_votes = hjson.loads(open(fullpath, 'r')) if os.path.exists(fullpath) else []

    for vote_id in votes_list:
        # Build url and fetch data for post
        full_url = urljoin(api_url, f"{vote_id}")
        data, response = get_data(full_url=full_url)

        # No data found for given post id
        if not data:
            logger.warning(f"Vote Id {vote_id} does not exist.")
            continue
        # Cache votes
        cached_votes.append(hjson.dumps(data))

    # Write to file.
    with open(fullpath, "w") as cache:
        hjson.dump(cached_votes, cache)
        logger.info(f"Cached {len(cached_votes)} votes into {fullpath}")

    return


def create_cache(root_dir, today=False, end_year=2009, start_year=2019):
    """
    Create user cache in root_dir from the start_date
    """
    # Get the remote site url
    api_url = "https://www.biostars.org/api/stats/date/"

    # Tuple with (year, month)
    months = [list(zip([year] * 12, range(1, 13))) for year in range(start_year, end_year + 1, -1)]
    # Flatten list
    months = list(chain(*months))
    # List of all days in given time period
    dates = [list(Calendar().itermonthdates(year=year, month=month)) for year, month in months]
    dates = list(chain(*dates))

    # Format the date to ensure its double digit.
    format_date = lambda digit: f"0{digit}" if len(f"{digit}") == 1 else f"{digit}"
    count = 0
    for date in dates:

        # Fetch the stats of the day
        month = format_date(date.month)
        day = format_date(date.day)
        daystr = f"{date.year}/{month}/{day}"
        full_url = urljoin(api_url, daystr)
        data, response = get_data(full_url=full_url)

        if not data:
            logger.error(f"No data for date={day}")
            continue

        # Create directory to store cache in
        prefix = f"{date.year}_{month}"
        os.makedirs(root_dir, exist_ok=True)

        # Get posts of the day
        cache_post(root_dir=root_dir, post_list=data.get("new_posts"), prefix=prefix)

        # Get users of the day
        cache_users(root_dir=root_dir, user_list=data.get("new_users"), prefix=prefix)

        # Get votes of the day
        cache_votes(root_dir=root_dir, votes_list=data.get("new_votes"), prefix=prefix)

        # TODO: take out after testing
        count += 1
        if count == 2:
            break
    return


def get_all_posts(post_files):
    """
    Return all posts in root dir as a dict keyed by post id.
    """
    all_posts = dict()

    for fname in post_files:
        posts = hjson.load(open(fname.path, "r"))
        for data in posts:
            post = hjson.loads(data)
            all_posts[str(post.get("id"))] = post

    return all_posts


def gen_post(post_id, all_posts):
    # Recursively create parent if it does not exist.
    post = all_posts[post_id]
    author_id = post.get("author_id")
    author = User.objects.filter(profile__uid=author_id).first()
    creation_date = post.get("creation_date")
    lastedit_date = post.get("lastedit_date")
    #has_accepted = post.get("has_accepted")
    lastedit_user_id = post.get("lastedit_user_id")
    lastedit_user = User.objects.filter(profile__uid=lastedit_user_id).first()

    parent_id = post.get("parent_id")
    root_id = post.get("root_id")
    parent = Post.objects.filter(uid=parent_id).first()
    root = Post.objects.filter(uid=root_id).first()
    title = post.get("title")
    post_type = post.get("type_id")
    view_count = post.get("view_count")
    tag_val = post.get("tag_val")
    html = post.get("xhtml")
    content = util.strip_tags(post.get("xhtml"))
    status = post.get("status_id")

    # Create root and parent if they do not exist.
    if not root_id and root_id != post_id:
        gen_post(post_id=parent_id, all_posts=all_posts)
    if not parent_id and parent_id != post_id:
        gen_post(post_id=parent_id, all_posts=all_posts)

    yield Post(content=content, uid=post_id, type=post_type, title=title, view_count=view_count,
               tag_val=tag_val, html=html, status=status, root=root, parent=parent, author=author,
               creation_date=creation_date, lastedit_date=lastedit_date, lastedit_user=lastedit_user)


def load_posts(post_files):

    all_posts = get_all_posts(post_files=post_files)

    def generate():
        # Iterate through all posts and create

        for post in all_posts:
            # Use function to recursively create parent and root posts that don't exists.
            gen_post(post_id=post, all_posts=all_posts)

    Post.objects.bulkcreate(objs=generate(), batch_size=50)
    return


def bulkload_votes(votes_file):

    def generate():

        votes = open(votes_file, "r").readlines()
        for data in votes:
            votes = hjson.loads(data)
            author_id = votes.get("author_id")
            post_id = votes.get("post_id")
            vote_type = votes.get("type_id")
            date = votes.get("date")
            uid = votes.get("uid")

            user = User.objects.filter(profile__uid=author_id)
            post = Post.objects.filter(uid=post_id)
            yield Vote(uid=uid, type=vote_type, date=date, author=user, post=post)

    # Bulk create votes
    Vote.objects.bulk_create(objs=generate(), batch_size=50)
    return


def bulkload_users(user_files):
    # TODO: ISSUE with bulk create, does not update profile.

    def generate():
        for fname in user_files:
            users = hjson.load(open(fname.path, "r"))
            for data in users:
                user = hjson.loads(data)
                #date_joined = user.get("date_joined")
                #last_login = user.get("last_login")
                name = user.get("name")
                uid = user.get("id")
                username = f"{name}{uid}"
                email = f"{name}{uid}@testmail.com"
                password = util.get_uuid()

                yield User(username=username, email=email, password=password)

    # Bulk create users.
    User.objects.bulk_create(objs=generate(), batch_size=50)
    return


def load_users(user_files):

    for fname in user_files:
        users = hjson.load(open(fname.path, "r"))
        for data in users:
            user = hjson.loads(data)
            date_joined = user.get("date_joined")
            last_login = user.get("last_login")
            name = user.get("name")
            uid = user.get("id")
            username = f"{name}{uid}"
            email = f"{name}{uid}@testmail.com"
            password = util.get_uuid()

            user = User.objects.filter(profile__uid=uid).first()
            if not user:
                user = User.objects.create(username=username, email=email, password=password)

            # Update the uid, last login and status
            Profile.objects.filter(user=user).update(uid=uid, date_joined=date_joined, last_login=last_login,
                                                      name=name)
            logger.info(f"Updated user={user.profile.uid}")

    return


def load_cache(root_dir):

    # Get the user, votes, and post files
    user_files = list(filter(lambda p: 'users' in p.name and p.is_file(), os.scandir(root_dir)))
    post_files = list(filter(lambda p: 'posts' in p.name and p.is_file(), os.scandir(root_dir)))
    votes_files = list(filter(lambda p: 'votes' in p.name and p.is_file(), os.scandir(root_dir)))

    # Load users first
    load_users(user_files=user_files)
    # Load posts next
    load_posts(post_files=post_files)
    # Bulk load votes last.
    bulkload_votes(votes_file=votes_files)

    return


class Command(BaseCommand):
    help = 'Cache posts and users from remote site.'

    def add_arguments(self, parser):
        parser.add_argument("--root", help="Root directory to store into or load from.")
        parser.add_argument("--load", action='store_true', help="Root directory to store into.")
        pass

    def handle(self, *args, **options):

        root = options["root"]
        load = options["load"]
        root = os.path.abspath(root)
        os.makedirs(root, exist_ok=True)

        if load:
            load_cache(root_dir=root)
        else:
            create_cache(root_dir=root)
        pass