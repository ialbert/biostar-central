from random import choice
from string import ascii_lowercase, digits

from django.contrib.auth import login, authenticate
from django.shortcuts import render, redirect
from .forms import SignUpForm
from django.contrib.auth import get_user_model


def already_used(username):

    usermodel = get_user_model()
    try:
        usermodel.objects.get(username=username)
        return True
    except usermodel.DoesNotExist:
        return False 


def hash_username(length=16):

    chars=ascii_lowercase+digits
    username = ''
    for token in range(length):
        username += choice(chars)

    if already_used(username):
        hash_username()
    else:
        return username 


def signup(request):

    if request.method == 'POST':
        form = SignUpForm(request.POST)

        # valid input at this point
        if form.is_valid():

            email = form.cleaned_data.get('email')
            raw_password = form.cleaned_data.get('password1')

            user = form.save(commit=False)
            user.username = hash_username()
            user.set_password(raw_password)
            user.email = email
            user.save()

            login(request, user)
            
            return redirect('/login')
    else:
        
        form = SignUpForm()
    return render(request, 'signup.html', {'form': form})


def index(request):
    return render(request, 'index.html')


#TEMPERORY ( not valid forms with csrf_tokens)
def login(request):
    return render (request, 'login.html')


def logout(request):
    return render(request, 'logout.html')



