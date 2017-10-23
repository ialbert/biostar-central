
import uuid
import logging
from django.contrib import messages

from ratelimit.decorators import ratelimit
from django.shortcuts import render, redirect
from django.contrib.auth import login, authenticate, logout

from .forms import SignUpForm, LoginForm, LogoutForm
from .models import User
from django.urls import reverse
from engine.const import *
from engine.views import breadcrumb_builder
from django.contrib import auth
from django.http import HttpResponseRedirect

logger = logging.getLogger('engine')


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]



def user_profile(request, id):

    user = User.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, USER_ICON], user=user)

    context = dict(user=user, steps=steps)

    return render(request, 'registration/user_profile.html', context)



@ratelimit(key='ip', rate='10/m', block=True, method=ratelimit.UNSAFE)
def user_signup(request):

    steps = breadcrumb_builder([HOME_ICON, SIGNUP_ICON])

    if request.method == 'POST':

        form = SignUpForm(request.POST)

        if form.is_valid():
            email = form.cleaned_data.get('email')
            password = form.cleaned_data.get('password1')
            name = email.split("@")[0]

            user = User.objects.create(username=get_uuid(), email=email,
                                       first_name=name)
            user.set_password(password)
            user.save()

            login(request, user)
            logger.info(f"Signed up and logged in user.id={user.id}, user.email={user.email}")
            messages.info(request, "Signup successful!")
            return redirect(reverse('login'))
    else:

        form = SignUpForm()
    context = dict(form=form, steps=steps)
    return render(request, 'registration/user_signup.html', context=context)


def user_logout(request):
    steps = breadcrumb_builder([HOME_ICON, LOGOUT_ICON])

    if request.method == "POST":

        form = LogoutForm(request.POST)

        if form.is_valid():
            auth.logout(request)
            messages.info(request, "You have been logged out")
            return redirect("/")


    form = LogoutForm()

    context = dict(steps=steps, form=form)

    return render(request, "registration/user_logout.html", context=context)


@ratelimit(key='ip', rate='10/m', block=True, method=ratelimit.UNSAFE)
def user_login(request):

    if request.method == "POST":
        form = LoginForm(data=request.POST)

        if form.is_valid():

            email = form.cleaned_data['email']
            password = form.cleaned_data['password']

            user = User.objects.filter(email=email).order_by('-id').first()

            if not user:
                form.add_error(None, "This email does not exist.")
            else:
                user = authenticate(username=user.username, password=password)
                login(request, user)
                logger.info(f"logged in user.id={user.id}, user.email={user.email}")
                messages.info(request, "Login successful!")

                return HttpResponseRedirect(request.GET.get('next', reverse("index")))

    else:
        initial = dict(nexturl=request.GET.get('next', '/'))
        form = LoginForm(initial=initial)

    steps = breadcrumb_builder([HOME_ICON, LOGIN_ICON])

    context = dict(form=form, steps=steps)
    return render(request, "registration/user_login.html", context=context)

