import json
import logging
from datetime import datetime, timedelta
from calendar import timegm
from os.path import join, normpath, isfile, exists
from os import makedirs

from django.http import HttpResponse
from django.contrib.sites.models import get_current_site
from django.conf import settings
from django.core.cache import get_cache

from ..apps.users.models import User
from ..apps.posts.models import Vote, Post, PostView


logger = logging.getLogger(__name__)


def json_response(f):
    """
    Converts any functions which returns a dictionary to a proper HttpResponse with json content.
    """
    # TODO: This decorator might be moved to util/json.py. We will see later on, when we
    # TODO: introduce django-rest-framework.
    def to_json(request, *args, **kwargs):
        """
        Creates the actual HttpResponse with json content.
        """
        data = f(request, *args, **kwargs)
        response = HttpResponse(json.dumps(data, sort_keys=True, indent=4),
                                content_type="application/json")
        if not data:
            response.status_code = 404
            response.reason_phrase = 'Not found'
        return response
    return to_json


@json_response
def traffic(request):
    """
    Traffic as post views in the last 60 min.
    """
    now = datetime.now()
    start = now-timedelta(minutes=60)
    try:
        post_views = PostView.objects.filter(date__gt=start).exclude(date__gt=now).distinct(
            'ip').count()
    except NotImplementedError:
        post_views = PostView.objects.filter(date__gt=start).exclude(date__gt=now).count()

    data = {
        'date': datetime_to_iso(now),
        'timestamp': datetime_to_unix(now),
        'post_views_last_60_min': post_views,
    }
    return data


@json_response
def user_details(request, id):
    """
    Details for a user.

    Parameters:
    id -- the id of the `User`.
    """
    try:
        user = User.objects.get(pk=id)
    except User.DoesNotExist:
        return {}

    days_ago = (datetime.now().date() - user.profile.date_joined.date()).days
    data = {
        'id': user.id,
        'name': user.name,
        'date_joined': datetime_to_iso(user.profile.date_joined),
        'last_login': datetime_to_iso(user.profile.last_login),
        'joined_days_ago': days_ago,
        'vote_count': Vote.objects.filter(author=user).count(),
    }
    return data


@json_response
def post_details(request, id):
    """
    Details for a post.

    Parameters:
    id -- the id of the `Post`.
    """
    try:
        post = Post.objects.get(pk=id)
    except Post.DoesNotExist:
        return {}

    data = {
        'id': post.id,
        'title': post.title,
        'type': post.get_type_display(),
        'type_id': post.type,
        'creation_date': datetime_to_iso(post.creation_date),
        'lastedit_date': datetime_to_iso(post.lastedit_date),
        'lastedit_user_id': post.lastedit_user.id,
        'author_id': post.author.id,
        'author': post.author.name,
        'status': post.get_status_display(),
        'status_id': post.status,
        'thread_score': post.thread_score,
        'rank': post.rank,
        'vote_count': post.vote_count,
        'view_count': post.view_count,
        'reply_count': post.reply_count,
        'comment_count': post.comment_count,
        'book_count': post.book_count,
        'subs_count': post.subs_count,
        'answer_count': post.root.reply_count,
        'has_accepted': post.has_accepted,
        'parent_id': post.parent.id,
        'root_id': post.root_id,
        'xhtml': post.html,
        'tag_val': post.tag_val,
        'url': 'http://{}{}'.format(get_current_site(request).domain, post.get_absolute_url()),
    }
    return data


@json_response
def vote_details(request, id):
    """
    Details for a vote.

    Parameters:
    id -- the id of the `Vote`.
    """
    try:
        vote = Vote.objects.get(pk=id)
    except Vote.DoesNotExist:
        return {}

    data = {
        'id': vote.id,
        'author_id': vote.author.id,
        'author': vote.author.name,
        'post_id': vote.post.id,
        'type': vote.get_type_display(),
        'type_id': vote.type,
        'date': datetime_to_iso(vote.date),
    }
    return data


@json_response
def daily_stats_on_day(request, day):
    """
    Statistics about this website for the given day.
    Day-0 is the day of the first post.

    Parameters:
    day -- a day, given as a number of days from day-0 (the day of the first post).
    """
    date = days_after_day_zero_to_datetime(day)

    # We don't provide stats for today or the future.
    if not date or date.date() >= datetime.today().date():
        return {}
    return compute_stats(date)


@json_response
def daily_stats_on_date(request, year, month, day):
    """
    Statistics about this website for the given date.

    Parameters:
    year -- Year, 4 digits.
    month -- Month, 2 digits.
    day -- Day, 2 digits.
    """
    date = datetime(int(year), int(month), int(day))
    # We don't provide stats for today or the future.
    if date.date() >= datetime.today().date():
        return {}
    return compute_stats(date)


## Statistics #####################################################################################

STATS_FOLDER = normpath(join(settings.EXPORT_DIR, 'stats'))


def compute_stats(date):
    """
    Statistics about this website for the given date.
    Statistics are stored to a json file for caching purpose.

    Parameters:
    date -- a `datetime`.
    """

    start = date.date()
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
        'date': datetime_to_iso(start),
        'timestamp': datetime_to_unix(start),
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

    if not settings.DEBUG:
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
    file_name = '{}-{}-{}.json'.format(date.year, date.month, date.day)
    return normpath(join(STATS_FOLDER, file_name))


## Date utils #####################################################################################

def datetime_to_iso(date):
    """
    Converts a datetime to the ISO8601 format, like: 2014-05-20T06:11:41.733900.

    Parameters:
    date -- a `datetime` instance.
    """
    if not isinstance(date, datetime):
        date = datetime.combine(date, datetime.min.time())
    return date.isoformat()


def datetime_to_unix(date):
    """
    Converts a datetime to a Unix timestamp , like: 1400566301.

    Parameters:
    date -- a `datetime` instance.
    """
    return timegm(date.timetuple())


def unix_to_datetime(timestamp):
    """
    Converts a Unix timestamp (like: 1400566301) to a datetime.

    Parameters:
    timestamp -- a Unix timestamp like: 1400566301.
    """
    return datetime.fromtimestamp(float(timestamp))


def days_after_day_zero_to_datetime(days):
    """
    Converts a date expressed as number of days after day-0 (the date of the first ever post) to
    `datetime`.

    Params:
    days -- number of days after day-0 (the date of the first post ever).
    """
    cache = get_cache('default')
    day_zero = cache.get('day_zero')

    if not day_zero:
        first_post = Post.objects.order_by('creation_date').only('creation_date')
        if not first_post:
            return False
        day_zero = first_post[0].creation_date
        cache.set('day_zero', day_zero, 60*60*24*7)  # Cache valid for a week.

    return day_zero + timedelta(days=int(days))