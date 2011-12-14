"""
Too many viewa in the main views.py

Started refactoring some here, this will eventually store all form based
actions whereas the main views.py will contain url based actions.
"""
from datetime import datetime, timedelta
from main.server import html, models, auth
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
from whoosh import index
from whoosh.qparser import QueryParser

class UserForm(forms.Form):
    "A form representing a new question"
    display_name = forms.CharField(max_length=30,  initial="", widget=forms.TextInput(attrs={'size':'30'}))   
    email        = forms.CharField(max_length=50,  initial="", widget=forms.TextInput(attrs={'size':'50'}))
    location     = forms.CharField(max_length=50,  required=False, initial="", widget=forms.TextInput(attrs={'size':'50'}))
    website      = forms.CharField(max_length=80,  required=False, initial="", widget=forms.TextInput(attrs={'size':'50'}))
    about_me     = forms.CharField(max_length=500, required=False, initial="", widget=forms.Textarea (attrs=dict(cols='50', rows=6)))

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
        

def user_edit(request, uid):
    "User's profile page"
    
    target = models.User.objects.select_related('profile').get(id=uid)
    
    allow = auth.authorize_user_edit(target=target, user=request.user, strict=False)
    if not allow:
        messages.error(request, "unable to edit this user")
        return html.redirect("/user/show/%s/" % uid)
        
    if request.method == 'GET':
        initial = dict(
            display_name = target.profile.display_name,
            email      = target.email or '',
            location   = target.profile.location or '',
            website    = target.profile.website or '',
            about_me   = target.profile.about_me or ''
        )
        form = UserForm(initial)
        return html.template(request, name='user.edit.html', user=target, form=form)
    elif request.method == 'POST':
        
        form = UserForm(request.POST)
        if not form.is_valid():
            return html.template(request, name='user.edit.html', user=target, form=form)
        else:
            for field in "display_name about_me website location".split():
                setattr(target.profile, field, form.cleaned_data[field])
            target.email = form.cleaned_data['email']
            target.profile.save()
            target.save()
            return html.redirect("/user/show/%s/" % target.id)

def about(request):
    "Renders the about page"

    post_count     = models.Post.objects.filter(status=POST_OPEN).count()
    question_count = models.Post.objects.filter(status=POST_OPEN, type=POST_QUESTION).count()
    answer_count   = models.Post.objects.filter(status=POST_OPEN, type=POST_ANSWER).count()
    comment_count  = models.Post.objects.filter(status=POST_OPEN, type=POST_COMMENT).count()
    user_count = models.User.objects.filter(profile__status=USER_ACTIVE).count()

    mods = models.User.objects.filter(profile__type=USER_MODERATOR).select_related("profile").all()[:100]
    admins = models.User.objects.filter(profile__type=USER_ADMIN).select_related("profile").all()[:100]

    params = html.Params(post_count=post_count, user_count=user_count, question_count=question_count, 
        answer_count=answer_count, comment_count=comment_count, admins=admins, mods=mods)

    return html.template(request, name='about.html', params=params)
   
def search(text):
    text = text.strip()[:200]
    if not text:
        return
    ix   = index.open_dir(settings.WHOOSH_INDEX)
    searcher = ix.searcher()
    parser   = QueryParser("content", ix.schema)
    query    = parser.parse(text)
    results  = searcher.search(query, limit=200)
    return [ hit['pid'] for hit in results ] 

def modlog_list(request):
    "Lists moderator actions"
    mods = models.Note.objects.filter(type=NOTE_MODERATOR).select_related('sender', 'target', 'post', 'sender_profile').order_by('-date')
    page = get_page(request, mods)
    return html.template(request, name='mod.log.list.html', page=page)

def badge_show(request, bid):
    "Shows users that have earned a certain badge"
    page = None
    badge  = models.Badge.objects.get(id=bid)
    awards = models.Award.objects.filter(badge=badge).select_related('user', 'user_profile')
    page  = get_page(request, awards, per_page=18)
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

@login_required(redirect_field_name='/openid/login/')
def destroy_post(request, pid):
    "Destroys a post"
    if request.method != 'POST':
        return html.json_response({'status':'error', 'msg':'Only POST requests are allowed'})        
   
    moderator = request.user
    post = models.Post.objects.get(id=pid)

    # for now only comments may be destroyed
    assert post.type == POST_COMMENT
    if not auth.authorize_post_edit(user=moderator, post=post, strict=False):
        return html.json_response({'status':'error', 'msg':'You do not have permission to delete this post.'})        
    
    post.delete()
    return html.json_response({'status':'success', 'msg':'post deleted'})
