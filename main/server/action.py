"""
Too many viewa in the main views.py

Started refactoring some here, this will eventually store all form based
actions whereas the main views.py will contain url based actions.
"""
from main.server import html, models
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
    first_name = forms.CharField(max_length=50,  initial="A")
    last_name  = forms.CharField(max_length=50,  initial="B")
    about_me   = forms.CharField(max_length=50,  initial="B", widget=forms.Textarea)
    
def user_edit(request, uid):
    "User's profile page"
    user = models.User.objects.select_related('profile').get(id=uid)
    
    if request.method == 'GET':
        initial = dict(
            first_name = user.first_name,
            last_name  = user.last_name,
            about_me   = user.profile.about_me
        )
        form = UserForm(initial)
    
        return html.template(request, name='user.edit.html', user=user, form=form)
      
