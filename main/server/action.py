"""
Too many viewa in the main views.py

Started refactoring some here, this will eventually store all form based
actions whereas the main views.py will contain url based actions.
"""
import os, sys, traceback

from datetime import datetime, timedelta
from main.server import html, models, auth, notegen
from main.server.html import get_page
from main.server.const import *

from django import forms
from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.contrib.auth import authenticate, login, logout
from django.conf import settings
from django.http import HttpResponse
from django.db.models import Q
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.core.mail import send_mail

from whoosh import index
from whoosh.qparser import QueryParser

# activate logging
import logging
logger = logging.getLogger(__name__)

class UserForm(forms.Form):
    "A form representing a new question"
    display_name = forms.CharField(max_length=30,  initial="", widget=forms.TextInput(attrs={'size':'30'}))   
    email        = forms.CharField(max_length=100,  initial="", widget=forms.TextInput(attrs={'size':'50'}))
    location     = forms.CharField(max_length=100,  required=False, initial="", widget=forms.TextInput(attrs={'size':'50'}))
    website      = forms.CharField(max_length=250,  required=False, initial="", widget=forms.TextInput(attrs={'size':'50'}))
    my_tags      = forms.CharField(max_length=250,  required=False, initial="", widget=forms.TextInput(attrs={'size':'50'}))
    about_me     = forms.CharField(max_length=2500, required=False, initial="", widget=forms.Textarea (attrs={'class':'span6'}))
    scholar      = forms.CharField(max_length=50,  required=False, initial="", widget=forms.TextInput(attrs={'size':'30'}))

LAST_CLEANUP = datetime.now()
def cleanup(request):
    "A call to this handler will attempt a database cleanup"
    global LAST_CLEANUP
    now  = datetime.now()
    diff = (now - LAST_CLEANUP).seconds
    if diff > 300: # five minutes        
        LAST_CLEANUP = now
        # get rid of unused tags
        models.Tag.objects.filter(count=0).delete()

@login_required(redirect_field_name='/openid/login/')
def post_moderate(request, pid, status):
    "General moderation function"
    user = request.user
    post = models.Post.objects.get(id=pid)
     
    # remap the status to valid
    status = dict(close=POST_CLOSED, open=POST_OPEN, delete=POST_DELETED).get(status)
    if not status:
        messages.error('Invalid post moderation action')
        return html.redirect( post.get_absolute_url() )    
    
    url = models.post_moderate(request=request, user=user, post=post, status=status)
    return html.redirect( url )    

@login_required(redirect_field_name='/openid/login/')
def user_moderate(request, uid, status):
    "General moderation function"
    user   = request.user
    target = models.User.objects.get(id=uid)
    url    = target.profile.get_absolute_url()

    # remap the status to valid
    status = dict(suspend=USER_SUSPENDED, reinstate=USER_ACTIVE).get(status)
    if not status:
        messages.error('Invalid user moderation action')
        return html.redirect( url )    
    
    flag, msg = models.user_moderate(user=user, target=target, status=status)
    func = messages.info if flag else messages.error
    func(request, msg)
    return html.redirect( url )    
    

@login_required(redirect_field_name='/openid/login/')
def user_edit(request, uid):
    "User's profile page"
    
    target = models.User.objects.select_related('profile').get(id=uid)
    
    allow = auth.authorize_user_edit(target=target, user=request.user, strict=False)
    if not allow:
        messages.error(request, "unable to edit this user")
        return html.redirect(target.profile.get_absolute_url() )
    
    # valid incoming fields
    fields = "display_name about_me website location my_tags scholar".split()
        
    if request.method == 'GET':
        initial = dict(email=target.email)
        for field in fields:
            initial[field] = getattr(target.profile, field) or ''                
        form = UserForm(initial)
        return html.template(request, name='user.edit.html', target=target, form=form)
    elif request.method == 'POST':
        
        form = UserForm(request.POST)
        if not form.is_valid():
            return html.template(request, name='user.edit.html', target=target, form=form)
        else:
            for field in fields:
                setattr(target.profile, field, form.cleaned_data[field])
            # check the new email
            new_email = form.cleaned_data['email'].strip()
            if new_email != target.email and models.User.objects.filter(email=new_email):
                # cannot set your email to an existing other user 'semail
                messages.error(request, "This email is aready taken - please merge the accounts!")
            else:
                target.email = new_email
            
            target.profile.save()
            target.save()
            
            url = reverse('main.server.views.user_profile', kwargs=dict(uid=target.id))
            return html.redirect(url)

@login_required(redirect_field_name='/openid/login/')
def post_reparent(request, pid, rid=0):
    "Reparent a post"
    
    post = models.Post.objects.get(id=pid)
    root = post.root
    parent = post.parent

    allow = auth.authorize_post_edit(post=post, user=request.user, strict=False)
    
    if not allow:
        messages.error(request, "Reparent access denied")
        return html.redirect(post.get_absolute_url())

    if post.type in POST_TOPLEVEL or post == post.root:
        messages.error(request, "Cannot reparent a toplevel post")
        return html.redirect(post.get_absolute_url())

    # these are the valid targets
    targets = models.Post.objects.filter(root=root).select_related('author', 'author__profile').exclude(id__in=(post.id, parent.id))

    target = request.REQUEST.get('target')
    if target:
        target =  models.Post.objects.get(id=target)
        
        if target not in targets:
            messages.error(request, "Invalid reparent %s -> %s" % (post.id, target.id) )
            return html.redirect(post.get_absolute_url())
        
        # comment to comment reparent is not yet supported
        if target.type == POST_COMMENT and post.type == POST_COMMENT:
            messages.error(request, "Comment to comment reparent %s -> %s not implemented" % (post.id, target.id) )
            return html.redirect(post.get_absolute_url())

        # perfomr the reparent
        post.parent = target
        question = (target.type == POST_QUESTION)
        post.type = POST_ANSWER if question else POST_COMMENT
        post.save()

        # valid target to be applied
        messages.info(request, "Reparenting %s to %s" % (post.id, target.id))
        return html.redirect(post.get_absolute_url())
        
    return html.template(request, name='post.reparent.html', post=post, targets=targets)

def badge_show(request, bid):
    "Shows users that have earned a certain badge"
    page = None
    badge  = models.Badge.objects.get(id=bid)
    awards = models.Award.objects.filter(badge=badge).select_related('user', 'user_profile').order_by("-date")
    page  = get_page(request, awards, per_page=24)
    return html.template(request, name='badge.show.html', page=page, badge=badge)
 
def note_clear(request, uid):
    "Clears all notifications of a user"
    user = models.User.objects.get(pk=uid)
    # you may only delete your own messages
    if user == request.user:
        messages.info(request, "All messages have been deleted")
        models.Note.objects.filter(target=user).all().delete()
    else:
        messages.warning(request, "You may only delete your own messages")
    return html.redirect("/user/show/%s/" % user.id)

ACCOUNT_MERGE_EMAIL = """

Account merge request by: %(request_name)s (%(request_id)s)

Master account: %(master_email)s, %(master_name)s
Master User: http://%(domain)s/user/profile/%(master_id)s/

Remove account: %(remove_email)s, %(remove_name)s
Remove User: http://%(domain)s/user/profile/%(remove_id)s/

Look at both user accounts to verify request.

To apply the merge click here:

http://%(domain)s/approve_merge/%(master_id)s/%(remove_id)s/

To ignore the request simply ignore this email.

"""

ACCOUNT_APPROVAL_EMAIL = """
Hello,

The requested BioStar account merge has been completed. 

Profile url: http://%(domain)s%(profile_url)s

cheers,

the BioStar Team
"""

@login_required(redirect_field_name='/openid/login/')
def request_merge(request):
    "Generates an account merge request"
    
    class MergeForm(forms.Form):
        "A form representing a new question"
        master_id = forms.CharField(max_length=5,  initial="", widget=forms.TextInput(attrs={'size':'5'}))   
        remove_id = forms.CharField(max_length=5,  initial="", widget=forms.TextInput(attrs={'size':'5'}))
    
    user = request.user
    
    if request.method == 'POST':   
        form = MergeForm(request.POST)
        if form.is_valid():
            try:
                data   = form.cleaned_data
                master = models.User.objects.get(id=data['master_id'])
                remove = models.User.objects.get(id=data['remove_id'])
                fill = dict(
                    domain=settings.SITE_DOMAIN, master_id=master.id, remove_id=remove.id,
                    master_name = master.profile.display_name, remove_name=remove.profile.display_name,
                    master_email = master.email, remove_email=remove.email,
                    request_id = request.user.id, request_name = user.profile.display_name,
                )
                body = ACCOUNT_MERGE_EMAIL % fill
                logger.info('sending email to %s' % settings.SERVER_EMAIL)
                send_mail('BioStar: account merge request', body, settings.DEFAULT_FROM_EMAIL, [ settings.SERVER_EMAIL ], fail_silently=False)
                messages.info(request, "Your request for account merge has been submitted for review.")
                return html.redirect( user.profile.get_absolute_url() )
            except Exception, exc:
                messages.error(request, 'Submission error %s' % exc)
    else: 
        form = MergeForm()
    params = html.Params(nav='')
    return html.template(request, name='pages/merge.html', params=params, form=form)

from django_openid_auth.models import UserOpenID

def migrate(master, remove):
    "Migrates user data"
    UserOpenID.objects.filter(user=remove).update(user=master)
    models.Vote.objects.filter(author=remove).update(author=master)
    models.Note.objects.filter(sender=remove).update(sender=master)
    models.Note.objects.filter(target=remove).update(target=master)
    models.Post.objects.filter(author=remove).update(author=master)
    models.Post.objects.filter(lastedit_user=remove).update(lastedit_user=master)
    models.PostRevision.objects.filter(author=remove).update(author=master)
    models.Award.objects.filter(user=remove).update(user=master)
    master.profile.score += remove.profile.score
    master.profile.save()
            
@login_required(redirect_field_name='/openid/login/')
def approve_merge(request, master_id, remove_id):
    "Approves an account merge request"
    user = request.user
    if not user.profile.is_admin:     
        messages.error(request, 'Error: approving user not an administrator!')
        return html.redirect("/")
   
    try:
        master = models.User.objects.get(id=master_id)
        remove = models.User.objects.get(id=remove_id)
        with transaction.commit_on_success():
            migrate(master, remove)
        remove.delete()
        fill = dict(
            domain=settings.SITE_DOMAIN, profile_url=master.profile.get_absolute_url()
        )
        body = ACCOUNT_APPROVAL_EMAIL % fill
        send_mail('BioStar: account merge complete', body, settings.DEFAULT_FROM_EMAIL, [ master.email ], fail_silently=False)
    except Exception, exc:
        messages.error(request, 'Merge error: %s' % exc)
        return html.redirect("/")
    
    messages.info(request, 'Merge completed')
    return html.redirect("/")

def test_login(request, uid, token):
    "This will allow test logins. Don't turn it on during production!"
    from django.contrib.auth import authenticate, login
    
    allow = (token == settings.SELENIUM_TEST_LOGIN_TOKEN)
    if settings.DEBUG and settings.SELENIUM_TEST_LOGIN_TOKEN and allow:
        user = models.User.objects.get(id=uid)
        password = models.make_uuid()
        user.set_password(password)
        user.save()
        user = authenticate(username=user.username, password=password)
        login(request=request, user=user)
        messages.info(request, "Test login complete.")
    else:
        messages.error(request, "Test login failed.")
        
    return html.redirect("/")   

def url500(request):
    "Custom error handler"
    
    type, value, tb = sys.exc_info()
    trace = traceback.format_exc()
    trace = "\n<trace>\n%s</trace>" % trace
    logger.error(trace)
    
    return html.template(request, name='500.html', path=request.path, value=value)
    
#
# this is only used to map redirects from the old site
#
POST_REMAP_FILE = '%s/db/post-remap.txt' % settings.HOME_DIR
if os.path.isfile(POST_REMAP_FILE):
    REMAP = dict( [line.split() for line in file(POST_REMAP_FILE)] )
else:
    REMAP = {}
    
def redirect_post(request, pid):
    try:
        nid = REMAP[pid]
        post = models.Post.objects.get(id=nid)
        return html.redirect(post.get_absolute_url(), permanent=True)   
    except Exception, exc:
        messages.error(request, "Unable to redirect: %s" % exc)
        return html.redirect("/")
        
def redirect_tag(request, tag):
    try:
        return html.redirect("/show/tag/%s/" % tag, permanent=True)   
    except Exception, exc:
        messages.error(request, "Unable to redirect: %s" % exc)
        return html.redirect("/")   