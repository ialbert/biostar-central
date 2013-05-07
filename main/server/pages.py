"""
semi-static pages 
"""
import urllib
from django.conf import settings
from main.server import html, models, formdef
from main.server.const import *
from django.db import connection
from django.contrib.sites.models import Site
from django.contrib import messages
from django.core.mail import send_mail
from datetime import datetime, timedelta
from django.core.urlresolvers import reverse


def about(request):
    "Renders the about page"

    recently = datetime.now() - timedelta(minutes=60)
    try:
        visitors = models.PostView.objects.filter(date__gt=recently).distinct('ip').count()
    except NotImplementedError, exc:
        visitors = models.PostView.objects.filter(date__gt=recently).count()

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
        answer_count=answer_count, comment_count=comment_count, admins=admins, mods=mods, navloc=navloc, managers=managers, visitors=visitors)
    
    return html.template(request, name='pages/about.html', params=params)
  
def rss(request):
    "Renders the rss feed page"
    user = request.user
    params = html.Params(nav='rss')
    return html.template(request, name='pages/rss.html', params=params, user=user)

def request_info(request, pid):
    "Requests information from a source"
    user = request.user
    post = models.Post.objects.get(id=pid)
    
    params = html.Params(site_domain = settings.SITE_DOMAIN, user=user, post=post)
    params.subject = "Your expert advice is needed at Biostar"
    
    if user.is_authenticated():
        params.display_name, score = user.profile.display_name, user.profile.score
    else:
        params.display_name, score = "Anonymous", 0
    
    params.body = html.fill(name='pages/request-info.txt', params=params)
     
    LIMIT = 5
    disabled = score < LIMIT
    
    if disabled:
        messages.error(request, "Note: users with fewer than %s reputation points may not send messages via Biostar. You have %s points" % (LIMIT, score))
    elif 'submit' in request.POST:
        form = formdef.RequestInfo(request.POST) 
        if form.is_valid():
            send_mail(subject=params.subject, from_email=settings.DEFAULT_FROM_EMAIL, message=params.body, recipient_list=[settings.DEFAULT_FROM_EMAIL, form.cleaned_data['email']], fail_silently=False)
            messages.info(request, "Your message has been sent.")
            return html.redirect( post.get_absolute_url() )
        else:
           messages.error(request, "%s" % form.errors) 
        
    return html.template(request, name='pages/request-info.html', params=params)
    
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


def help_external(request, word='main'):
    "Renders a test page"
    user = request.user

    # create a test key/value pair
    name = "TEST-KEY"
    pair = settings.EXTERNAL_AUTHENICATION.get(name)
    if pair:
        messages.info(request, "The <b>TEST KEY</b> is active!")
        key, patt = pair
    else:
        key, patt = "abcd", ""

    # prepare the data
    data = dict(display_name="Jane Doe",
                username="galaxy-jane2",
                email="jane.doe@gmail.com", tags="galaxy bwa", title="How do I run bwa?",
                tool_name="Uber Convert",
                tool_id=1,
                tool_version="0.0.1.alpha",
    )

    # encode the data and get the digest
    enc, digest = formdef.encode(data, key=key)

    # encode into url parameters
    store = dict(name=name, data=enc, digest=digest)

    # encoding the parameters into the url to be loaded
    params = urllib.urlencode(store.items())
    login_url = "/x/?%s" % params

    store['action'] = 'new'
    params = urllib.urlencode(store.items())
    post_url = "/x/?%s" % params

    # this is used inside the templates
    params = html.Params(post_url=post_url, login_url=login_url, key=key, data=data, enc=enc, digest=digest)
    return html.template(request, name='help/help.external.html', params=params, user=user)

def testpage(request):
    "Renders a test page"
    user = request.user
    params = html.Params()
    return html.template(request, name='pages/testpage.html', params=params, user=user)