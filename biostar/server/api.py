import json
from datetime import datetime, timedelta

from django.http import HttpResponse
from django.contrib.sites.models import get_current_site

from .util.dates import datetime_to_iso, datetime_to_unix
from .util.stats import count_traffic, compute_stats
from ..apps.users.models import User
from ..apps.posts.models import Vote, Post


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
    data = {
        'date': datetime_to_iso(now),
        'timestamp': datetime_to_unix(now),
        'post_views_last_60_min': count_traffic(now-timedelta(minutes=60), now),
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
        'author_id ': vote.author.id,
        'author ': vote.author.name,
        'post_id': vote.post.id,
        'type': vote.get_type_display(),
        'type_id': vote.type,
        'date': datetime_to_iso(vote.date),
    }
    return data


@json_response
def daily_stats(request, days_ago=1):
    """
    Statistics about this website for the given day.

    Parameters:
    days_ago -- a day, given as a number of days ago.
    """
    days_ago = int(days_ago)

    if days_ago < 1:
        return {}

    return compute_stats(days_ago)