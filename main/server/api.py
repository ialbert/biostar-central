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

        resp = HttpResponse(simplejson.dumps(data, sort_keys=True, indent=4), mimetype="application/json")
        return resp

def date_interval(days):
    """Returns the date limits that frame the date "days" ago"""
    days = int(days)
    now  = datetime.now().date()
    start = now - timedelta(days=days)
    end = start + timedelta(days=1)
    return start, end

def user_wrap(user):
    now  = datetime.now().date()
    ago  = (now - user.date_joined.date()).days
    data =  dict(id=user.id, name=user.profile.display_name,
        date_joined = timefmt(user.date_joined),
        last_visited= timefmt(user.profile.last_visited),
        joined_days_ago = ago,
    )
    return data

def post_wrap(post):
    data =  dict(id=post.id, title=post.title,
        type=post.get_type_display(),
        type_id=post.type,
        creation_date = timefmt(post.creation_date),
        lastedit_date= timefmt(post.lastedit_date),
        author_id = post.author.id,
        author = post.author.profile.display_name,
        score = post.score,
        answer_count = post.answer_count,
        rank = post.rank,
        parent_id = post.parent.id,
        root_id=post.root_id,
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

def join_path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

API_DIR = join_path(settings.EXPORT_DIR, 'api')
if not os.path.isdir(API_DIR):
    print '*** generating API_DIR=%s' % API_DIR
    os.mkdir(API_DIR)

def stats(request, days=1):
    days = int(days)

    global API_DIR

    if days < 1:
        return HttpResponse("", mimetype="application/json")

    start, end = date_interval(days=days)

    fname = time.strftime("%Y-%m-%d.json", start.timetuple())
    fpath = join_path(API_DIR, fname)

    if not os.path.isfile(fpath):
        data = make_stats(days)
        json = simplejson.dumps(data, sort_keys=True, indent=4)
        print '*** writing api file=%s' % fpath
        fp = file(fpath, 'wt')
        fp.write(json)
        fp.close()

    # return the precomputed data
    json = file(fpath).read()
    return HttpResponse(json, mimetype="application/json")

def make_stats(days=0):
    "This returns the json data about biostar"

    start, end = date_interval(days=days)


    now  = datetime.now().date()
    ago  = (now - start).days

    query = models.Post.objects.filter

    questions = query(type=POST_QUESTION, creation_date__lt=end).count()
    answers   = query(type=POST_ANSWER, creation_date__lt=end).count()
    toplevel  = query(type__in=POST_TOPLEVEL, creation_date__lt=end).exclude(type=POST_BLOG).count()
    comments  = query(type=POST_COMMENT, creation_date__lt=end).count()
    votes     = models.Vote.objects.filter(date__lt=end).count()
    users     = models.User.objects.filter(date_joined__lt=end).count()

    # out the x here so that these are displayed at the end
    x_new_users = models.User.objects.filter(date_joined__gte=start, date_joined__lt=end)
    x_new_users = map(user_wrap, x_new_users)

    x_new_posts = models.Post.objects.filter(creation_date__gte=start, creation_date__lt=end)
    x_new_posts = map(post_wrap, x_new_posts)

    x_new_votes = models.Vote.objects.filter(date__gte=start, date__lt=end)
    x_new_votes = map(vote_wrap, x_new_votes)

    data = {
        'days_ago': ago,
        'date': timefmt(start),
        'timestamp': time.mktime(end.timetuple()),
        'questions': questions,
        'answers': answers,
        'toplevel': toplevel,
        'comments': comments,
        'votes': votes,
        'users': users,
        'x_new_users': x_new_users,
        'x_new_posts': x_new_posts,
        'x_new_votes': x_new_votes,
    }
    return data
