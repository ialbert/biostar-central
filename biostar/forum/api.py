
import json
import os
import logging
from os.path import join, normpath
from django.core.cache import cache
from django.conf import settings
from datetime import datetime, timedelta

from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt
from biostar.accounts.models import Profile, User
from . import util
from .models import Post, Vote, Subscription, PostView


logger = logging.getLogger("engine")


def api_error(msg="Api Error"):
    return {'error': msg}


def stat_file(date, data=None, load=False, dump=False):

    os.makedirs(settings.STATS_DIR, exist_ok=True)
    file_name = f'{date.year}-{date.month}-{date.day}.json'
    file_path = normpath(join(settings.STATS_DIR, file_name))

    def load_file():
        # This will be FileNotFoundError in Python3.
        if not os.path.isfile(file_path):
            raise IOError
        with open(file_path, 'r') as fin:
            return json.loads(fin.read())

    def dump_into_file():
        with open(file_path, 'w') as fout:
            fout.write(json.dumps(data))

    if load:
        return load_file()

    if dump:
        return dump_into_file()


def get_counts(end):
    questions = Post.objects.filter(type=Post.QUESTION, creation_date__lt=end).count()
    answers = Post.objects.filter(type=Post.ANSWER, creation_date__lt=end).count()
    toplevel = Post.objects.filter(type__in=Post.TOP_LEVEL, creation_date__lt=end).exclude(type=Post.BLOG).count()
    comments = Post.objects.filter(type=Post.COMMENT, creation_date__lt=end).count()
    votes = Vote.objects.filter(date__lt=end).count()
    users = User.objects.filter(profile__date_joined__lt=end).count()

    data = {
        'questions': questions,
        'answers': answers,
        'toplevel': toplevel,
        'comments': comments,
        'votes': votes,
        'users': users,
    }
    return data


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
        return stat_file(date=start, load=True)
    except Exception as exc:  # This will be FileNotFoundError in Python3.
        logger.info('No stats file for {}.'.format(start))

    new_users = Profile.objects.filter(date_joined__gte=start,
                                       date_joined__lt=end).values_list("uid", flat=True)
    new_posts = Post.objects.filter(creation_date__gte=start,
                                    creation_date__lt=end).values_list("uid", flat=True)
    new_votes = Vote.objects.filter(date__gte=start,
                                    date__lt=end).values_list("id", flat=True)

    data = {
        'date': util.datetime_to_iso(start),
        'timestamp': util.datetime_to_unix(start),
        'new_users': list(new_users),
        'new_posts': list(new_posts),
        'new_votes': list(new_votes),
    }

    data.update(get_counts(end=end))

    if not settings.DEBUG:
        stat_file(dump=True, date=start, data=data)

    return data


def json_response(f):
    """
    Converts any functions which returns a dictionary to a proper HttpResponse with json content.
    """
    def to_json(request, *args, **kwargs):
        """
        Creates the actual HttpResponse with json content.
        """
        try:
            data = f(request, *args, **kwargs)
        except Exception as exc:
            logger.error(exc)
            data = api_error(msg=f"Error: {exc}")

        payload = json.dumps(data, sort_keys=True, indent=4)
        response = HttpResponse(payload, content_type="application/json")
        if not data:
            response.status_code = 404
            response.reason_phrase = 'Not found'
        return response
    return to_json


@json_response
def daily_stats_on_day(request, day):
    """
    Statistics about this website for the given day.
    Day-0 is the day of the first post.

    Parameters:
    day -- a day, given as a number of days from day-0 (the day of the first post).
    """
    day_zero = cache.get('day_zero')

    first_post = Post.objects.order_by('creation_date').only('creation_date')

    if day_zero is None and not first_post:
        return False

    if day_zero is None:
        day_zero = first_post[0].creation_date
        cache.set('day_zero', day_zero, 60 * 60 * 24 * 7)  # Cache valid for a week.

    date = day_zero + timedelta(days=int(day))

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


@json_response
def traffic(request):
    """
    Traffic as post views in the last 60 min.
    """
    now = datetime.now()
    start = now - timedelta(minutes=60)
    start = start.now()
    post_views = PostView.objects.filter(date__gt=start).exclude(date__gt=now).distinct('ip').count()

    data = {
        'date': util.datetime_to_iso(now),
        'timestamp': util.datetime_to_unix(now),
        'post_views_last_60_min': post_views,
    }
    return data


@json_response
def api_tag(request, tag):
    """
    Return list of post uids that have a tag.
    """
    posts = Post.objects.filter(tags__name=tag.lower()).values_list('uid', flat=True)
    posts = list(posts)

    return posts


@json_response
def user_email(request, email):
    user = User.objects.filter(email=email.lower())
    if user.exists():
        return True

    return False


@json_response
def user_details(request, uid):
    """
    Details for a user.

    Parameters:
    id -- the uid of the `User`.
    """

    user = User.objects.filter(profile__uid=uid).first()
    if not user:
        return {}

    days_ago = (datetime.now().date() - user.profile.date_joined.date()).days
    data = {
        'uid': user.profile.uid,
        'name': user.profile.name,
        'date_joined': util.datetime_to_iso(user.profile.date_joined),
        'last_login': util.datetime_to_iso(user.profile.last_login),
        'joined_days_ago': days_ago,
        'vote_count': Vote.objects.filter(author=user).count(),
    }
    return data


@json_response
def post_details(request, uid):
    """
    Details for a post.

    Parameters:
    id -- the id of the `Post`.
    """

    post = Post.objects.filter(uid=uid).first()
    if not post:
        return {}
    return post.json_data()


@json_response
def watched_tags(request, email):
    """
    Show watched tags for a user, given API key.
    Parameters:
    uid -- the id of the `User`.
    """
    user = User.objects.filter(email=email.lower()).first()
    if user:
        data = {'watched_tags': user.profile.watched_tags}
    else:
        data = {}

    return data


@json_response
def vote_details(request, uid):
    """
    Details for a vote.

    Parameters:
    uid -- the id of the `Vote`.
    """
    vote = Vote.objects.filter(id=uid).first()

    if not vote:
        return {}

    data = {
        'id': vote.id,
        'author_uid': vote.author.profile.uid,
        'author': vote.author.profile.name,
        'post_uid': vote.post.uid,
        'type': vote.get_type_display(),
        'type_id': vote.type,
        'date': util.datetime_to_iso(vote.date),
    }
    return data


@csrf_exempt
@json_response
def tags_list(request):
    """
    Given a file of tags, return the post count for each.
    """
    # Get a file with all the tags.
    tags = request.FILES.get('tags')

    # How many months prior to look back
    months = request.POST.get('months', '6')

    try:
        months = int(months) if months.isalnum() else float(months)
    except Exception as exc:
        logger.error(exc)
        months = 6

    # Convert months to weeks
    weeks = months * 4

    delta = util.now() - timedelta(weeks=weeks)
    query = Post.objects.filter(lastedit_date__gt=delta)

    # Iterate over tags and collect counts.
    lines = tags.readlines() if tags else []

    data = {}

    for line in lines:
        tag = line.decode().lower().strip()
        if not tag:
            continue

        # Get posts with this tag, that have been answered
        posts = query.filter(tags__name=tag.lower(), is_toplevel=True)
        uids = posts.values_list('uid', flat=True)

        # Filter so the count is only 1 answer per questions
        answer_count = Post.objects.filter(uid__in=uids, answer_count__gte=1).count()
        comment_count = Post.objects.filter(uid__in=uids, comment_count__gte=1).count()
        total = len(uids)

        val = dict(total=total, answer_count=answer_count, comment_count=comment_count)
        data.setdefault(tag, {}).update(val)

    return data
