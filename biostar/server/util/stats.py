"""
Utility module for statistics about the website.
"""
from datetime import datetime, timedelta
import logging
import time
from os.path import join, normpath, isfile, exists
from os import makedirs
import json

from django.conf import settings

from biostar.apps.posts.models import PostView, Post, Vote
from biostar.apps.users.models import User
from .dates import datetime_to_iso, datetime_to_unix


logger = logging.getLogger(__name__)


STATS_FOLDER = normpath(join(settings.EXPORT_DIR, 'stats'))


def count_traffic(start, end):
    """
    Returns the traffic as post views in the time interval defined by `start` and `end`.

    Parameters:
    start -- a `datetime` instance.
    end -- a `datetime` instance.
    """
    try:
        return PostView.objects.filter(date__gt=start).exclude(date__gt=end).distinct('ip').count()
    except NotImplementedError:
        return PostView.objects.filter(date__gt=start).exclude(date__gt=end).count()


def compute_stats(days_ago=0):
    """
    Statistics about this website for the given day.
    Statistics are stored to a json file for caching purpose.

    Parameters:
    days_ago -- a day, given as a number of days ago.
    """

    now = datetime.now().date()
    start = now - timedelta(days=days_ago)
    end = start + timedelta(days=1)

    try:
        return load_stats_from_file(start)
    except IOError:  # This will be FileNotFoundError in Python3.
        logger.info('No stats file for {}.'.format(start))

    query = Post.objects.filter

    questions = query(type=Post.QUESTION, creation_date__lt=end).count()
    answers = query(type=Post.ANSWER, creation_date__lt=end).count()
    toplevel = query(type__in=Post.TOP_LEVEL, creation_date__lt=end).exclude(
        type=Post.BLOG).count()
    comments = query(type=Post.COMMENT, creation_date__lt=end).count()
    votes = Vote.objects.filter(date__lt=end).count()
    users = User.objects.filter(profile__date_joined__lt=end).count()

    new_users = User.objects.filter(profile__date_joined__gte=start, profile__date_joined__lt=end)
    new_users_ids = [user.id for user in new_users]

    new_posts = Post.objects.filter(creation_date__gte=start, creation_date__lt=end)
    new_posts_ids = [post.id for post in new_posts]

    new_votes = Vote.objects.filter(date__gte=start, date__lt=end)
    new_votes_ids = [vote.id for vote in new_votes]

    data = {
        'days_ago': days_ago,
        'date': datetime_to_iso(now),
        'timestamp': datetime_to_unix(end),
        'questions': questions,
        'answers': answers,
        'toplevel': toplevel,
        'comments': comments,
        'votes': votes,
        'users': users,
        'new_users': new_users_ids,
        'new_posts': new_posts_ids,
        'new_votes': new_votes_ids,
    }
    dump_stats_to_file(start, data)
    return data


def load_stats_from_file(date):
    """
    Load stats data from a stat file.

    Params:
    date -- a `datetime` instance.
    """
    file_path = _build_stats_file_path(date)

    if not isfile(file_path):
        raise IOError  # This will be FileNotFoundError in Python3.

    with open(file_path, 'r') as fin:
        return json.loads(fin.read())


def dump_stats_to_file(date, data):
    """
    Store stats data to a json file for caching purpose.

    Params:
    date -- a `datetime` instance.
    data -- stats data in a dictionary.
    """
    if not exists(STATS_FOLDER):
        # Ensure STATS_FOLDER exists.
        makedirs(STATS_FOLDER)

    file_path = _build_stats_file_path(date)

    with open(file_path, 'w') as fout:
        fout.write(json.dumps(data))


def _build_stats_file_path(date):
    """
    Build the path of a stat file from a date.

    Params:
    date -- a `datetime` instance.
    """
    file_name = time.strftime("%Y-%m-%d.json", date.timetuple())
    return normpath(join(STATS_FOLDER, file_name))