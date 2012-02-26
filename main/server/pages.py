"""
semi-static pages 
"""
from main.server import html, models
from main.server.const import *

def about(request):
    "Renders the about page"

    post_count     = models.Post.objects.filter(status=POST_OPEN).count()
    question_count = models.Post.objects.filter(status=POST_OPEN, type=POST_QUESTION).count()
    answer_count   = models.Post.objects.filter(status=POST_OPEN, type=POST_ANSWER).count()
    comment_count  = models.Post.objects.filter(status=POST_OPEN, type=POST_COMMENT).count()
    user_count = models.User.objects.filter(profile__status=USER_ACTIVE).count()

    mods = models.User.objects.filter(profile__type=USER_MODERATOR).select_related("profile").all()[:100]
    admins = models.User.objects.filter(profile__type=USER_ADMIN).select_related("profile").all()[:100]

    navloc = dict(about="active")
    params = html.Params(nav='about', post_count=post_count, user_count=user_count, question_count=question_count, 
        answer_count=answer_count, comment_count=comment_count, admins=admins, mods=mods, navloc=navloc)
    
    return html.template(request, name='about.html', params=params)
  
def rss(request):
    "Renders the rss feed page"
    params = html.Params(nav='rss')
    return html.template(request, name='rss.html', params=params)

def faq(request):
    "Renders the faq page"
    params = html.Params(nav='faq')
    return html.template(request, name='faq.html', params=params)

def beta(request):
    "Renders the beta test information page"
    params = html.Params(nav='')
    return html.template(request, name='beta.html', params=params)

  