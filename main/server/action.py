"""
Too many viewa in the main views.py

Started refactoring some here, this will eventually store all form based
actions whereas the main views.py will contain url based actions.
"""
from main.server import html, models
from main.server.html import get_page

from django import forms
from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.contrib.auth import authenticate, login, logout
from django.conf import settings
from django.http import HttpResponse
from django.db.models import Q

class UserForm(forms.Form):
    "A form representing a new question"
    display_name = forms.CharField(max_length=30,  initial="B", widget=forms.TextInput(attrs={'size':'30'}))     
    location     = forms.CharField(max_length=50,  initial="B", widget=forms.TextInput(attrs={'size':'50'}))
    website      = forms.CharField(max_length=80,  initial="B", widget=forms.TextInput(attrs={'size':'50'}))
    about_me     = forms.CharField(max_length=500, initial="C", widget=forms.Textarea (attrs=dict(cols='50', rows=6)))
    
def user_edit(request, uid):
    "User's profile page"
    user = models.User.objects.select_related('profile').get(id=uid)
    
    if request.method == 'GET':
        initial = dict(
            display_name = user.profile.display_name,
            location   = user.profile.location or 'not specified',
            website    = user.profile.website or 'http://www.biostars.org',
            about_me   = user.profile.about_me or 'not specified'
        )
        form = UserForm(initial)
        return html.template(request, name='user.edit.html', user=user, form=form)
    elif request.method == 'POST':
        
        form = UserForm(request.POST)
        if not form.is_valid():
            return html.template(request, name='user.edit.html', user=user, form=form)
        else:
            for field in "display_name about_me website location".split():
                setattr(user.profile, field, form.cleaned_data[field])
            user.profile.save()
            user.save()
            return html.redirect("/user/show/%s/" % user.id)
         
def modlog_list(request):
    "Lists moderator actions"
    mods = models.ModLog.objects.order_by('-date')
    page = get_page(request, mods)
    return html.template(request, name='modlog.list.html', page=page)
    
        