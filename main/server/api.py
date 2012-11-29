"""
Too many viewa in the main views.py

Started refactoring some here, this will eventually store all form based
actions whereas the main views.py will contain url based actions.
"""
import os, sys, traceback, time
from datetime import datetime, timedelta
from main.server import html, models
from main.server.const import *
from django.conf import settings
from django.http import HttpResponse
from django.db.models import Q
from django.utils import simplejson
from django.utils import simplejson as json

# earliest date is 1155 days ago

# activate logging
import logging
logger = logging.getLogger(__name__)


class json_wrapper(object):
    "Used as decorator to trap  errors in json calls"
    def __init__(self, f):
        self.f = f

    def __call__(self, *args, **kwds):
        try:
            now = datetime.now()
            data = self.f(*args, **kwds)
        except Exception,exc:
            traceback.print_exc()
            data = dict(error='%s' % exc)

        resp = HttpResponse(json.dumps(data, sort_keys=True, indent=4), mimetype="application/json")
        return resp

def date_interval(days):
    """Returns the date limits that frame the date "days" ago"""
    days = int(days)
    now  = datetime.now()
    past = now - timedelta(days=days)
    start = past.date()
    end = start + timedelta(days=1)
    return start, end

def user_wrap(user):
    data =  dict(id=user.id, name=user.profile.display_name,
        date_joined = timefmt(user.date_joined),
        last_visited= timefmt(user.profile.last_visited)
    )
    return data

def post_wrap(post):
    data =  dict(id=post.id, title=post.title,
        type=post.get_type_display(),
        creation_date = timefmt(post.creation_date),
        lastedit_date= timefmt(post.lastedit_date),
        author_id = post.author.id,
        author = post.author.profile.display_name,
        score = post.score,
        answer_count = post.answer_count,
        rank = post.rank,
        parent_id = post.parent.id,
    )
    return data

def vote_wrap(vote):
    data =  dict(
        post_id=vote.post.id,
        date=timefmt(vote.date),
        type=vote.get_type_display(),
    )
    return data

@json_wrapper
def user_info(request, uid):
    user = models.User.objects.get(pk=uid)
    data = user_wrap(user)
    data['vote_count'] = models.Vote.objects.filter(author=user).count()
    return data

@json_wrapper
def post_info(request, pid):
    post = models.Post.objects.get(pk=pid)
    data = post_wrap(post)
    data['xhtml'] = post.html
    return data

@json_wrapper
def new_posts(request, days):
    start, end = date_interval(days)
    posts = models.Post.objects.filter(creation_date__gte=start, creation_date__lt=end)
    post_count = models.Post.objects.filter(creation_date__lt=end).count()
    date = timefmt(start)

    new_posts = map(post_wrap, posts)
    data = dict(
        description="posts created on %s" % date,
        new_posts = new_posts,
        date=date,
        post_count = post_count,
    )
    return data

@json_wrapper
def new_votes(request, days):
    start, end = date_interval(days)
    votes = models.Vote.objects.filter(date__gte=start, date__lt=end)
    vote_count = models.Vote.objects.filter(date__lt=end).count()
    date = timefmt(start)

    new_votes = map(vote_wrap, votes)
    data = dict(
        description="votes created on %s" % date,
        new_votes = new_votes,
        date=date,
        vote_count = vote_count,
    )
    return data

@json_wrapper
def new_users(request, days=0):
    start, end = date_interval(days)
    users = models.User.objects.filter(date_joined__gte=start, date_joined__lt=end)
    user_count = models.User.objects.filter(date_joined__lt=end).count()
    new_users = map(user_wrap, users)
    date = timefmt(start)
    data = dict(new_users=new_users, user_count=user_count, date=date)
    return data

def timefmt(date):
    return time.strftime("%Y-%m-%d", date.timetuple())

def get_traffic(end, minutes=60):
    "Returns the traffic as a number"
    try:
        start = end - timedelta(minutes=minutes)
        traffic = models.PostView.objects.filter(date__gt=start).exclude(date__gt=end).distinct('ip').count()
    except NotImplementedError, exc:
        traffic = models.PostView.objects.filter(date__gt=start).exclude(date__gt=end).count()
    return traffic

@json_wrapper
def traffic(request):
    now = datetime.now()
    minutes = 60;
    data = {
        'date': now.ctime(),
        'timestamp': time.mktime(now.timetuple()),
        'traffic': get_traffic(now, minutes=minutes),
    }
    return data

@json_wrapper
def stats(request, days=0):
    "This return a json data about biostar"

    now = datetime.now()
    start, end = date_interval(days=days)

    query = models.Post.objects.filter
    minutes = 60;
    data = {
        'date': end.ctime(),
        'timestamp': time.mktime(end.timetuple()),
        'questions': query(type=POST_QUESTION, creation_date__lt=end).count(),
        'answers': query(type=POST_ANSWER, creation_date__lt=end).count(),
        'toplevel': query(type__in=POST_TOPLEVEL, creation_date__lt=end).exclude(type=POST_BLOG).count(),
        'comments': query(type=POST_COMMENT, creation_date__lt=end).count(),
        'votes':  models.Vote.objects.filter(date__lt=end).count(),
        'users': models.User.objects.filter(date_joined__lt=end).count(),
    }
    return data
