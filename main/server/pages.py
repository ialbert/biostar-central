"""
semi-static pages 
"""
from django.conf import settings
from main.server import html, models
from main.server.const import *
from django.db import connection
from django.contrib.sites.models import Site

def about(request):
    "Renders the about page"

    post_count     = models.Post.objects.filter(status=POST_OPEN).count()
    question_count = models.Post.objects.filter(status=POST_OPEN, type=POST_QUESTION).count()
    answer_count   = models.Post.objects.filter(status=POST_OPEN, type=POST_ANSWER).count()
    comment_count  = models.Post.objects.filter(status=POST_OPEN, type=POST_COMMENT).count()
    user_count = models.User.objects.filter(profile__status=USER_ACTIVE).count()

    mods = models.User.objects.filter(profile__type=USER_MODERATOR).select_related("profile").order_by('-profile__score').all()
    admins = models.User.objects.filter(profile__type=USER_ADMIN).select_related("profile").order_by('-profile__score').all()
    managers = models.User.objects.filter(email=settings.ADMINS[0][1]).select_related("profile").order_by('-profile__score').all()
    navloc = dict(about="active")
    params = html.Params(nav='about', post_count=post_count, user_count=user_count, question_count=question_count, 
        answer_count=answer_count, comment_count=comment_count, admins=admins, mods=mods, navloc=navloc, managers=managers)
    
    return html.template(request, name='pages/about.html', params=params)
  
def rss(request):
    "Renders the rss feed page"
    user = request.user
    params = html.Params(nav='rss')
    return html.template(request, name='pages/rss.html', params=params, user=user)

def google(request):
    "Renders the rss feed page"
    user = request.user
    params = html.Params(nav='google')
    return html.template(request, name='pages/google.html', params=params, user=user)
    
def faq(request):
    "Renders the faq page"
    best = models.User.objects.all().select_related("profile").order_by('-profile__score')[:3]
    params = html.Params(nav='faq', best=best)
    return html.template(request, name='pages/faq.html', params=params)

def beta(request):
    "Renders the beta test information page"
    params = html.Params(nav='')
    return html.template(request, name='pages/beta.html', params=params)

def testpage(request):
    "Renders a test page"
    user = request.user
    params = html.Params(nav='rss')

    posts = models.Post.objects
    posts = posts.exclude(type=POST_BLOG).select_related('author', 'author__profile')
    posts = posts.order_by('-rank')[:10]
    posts = list(posts)


    rows = connection.queries

    return html.template(request, name='pages/testpage.html', params=params, user=user, rows=rows)
