
import uuid
import logging

from django.utils import http
from django.contrib import messages
from ratelimit.decorators import ratelimit
from django.views.decorators import csrf, cache
from django.shortcuts import render, redirect
from django.contrib.auth.models import User
from django.contrib.auth import views as auth_views
from django.urls import reverse


from .forms import SignUpForm, LoginForm, LogoutForm
from biostar.engine.const import *
from biostar.engine.views import breadcrumb_builder
from django.contrib import auth

logger = logging.getLogger('engine')

def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def accounts(request):

    return redirect("/")


def user_profile(request, id):

    user = User.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, USER_ICON], user=user)

    context = dict(user=user, steps=steps)

    return render(request, 'accounts/profile.html', context)



@ratelimit(key='ip', rate='10/m', block=True, method=ratelimit.UNSAFE)
@csrf.csrf_protect
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

            auth.login(request, user)
            logger.info(f"Signed up and logged in user.id={user.id}, user.email={user.email}")
            messages.info(request, "Signup successful!")
            return redirect(reverse('login'))
    else:

        form = SignUpForm()
    context = dict(form=form, steps=steps)
    return render(request, 'accounts/signup.html', context=context)


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

    return render(request, "accounts/logout.html", context=context)


@ratelimit(key='ip', rate='10/m', block=True, method=ratelimit.UNSAFE)
@csrf.csrf_protect
@cache.never_cache
def user_login(request):

    steps = breadcrumb_builder([HOME_ICON, LOGIN_ICON])

    if request.method == "POST":
        auth.logout(request)
        form = LoginForm(data=request.POST)
        next = request.POST.get('next', request.GET.get("next", "/"))
        if not http.is_safe_url(next, request.get_host()):
            next = "/"

        if form.is_valid():

            email = form.cleaned_data['email']
            password = form.cleaned_data['password']
            # Due to an early bug emails may not be unique. Last subscription wins.
            user = User.objects.filter(email__iexact=email).order_by('-id').first()

            if not user:
                form.add_error(None, "This email does not exist.")
                context = dict(form=form, steps=steps)
                return render(request, "accounts/login.html", context=context)

            user = auth.authenticate(username=user.username, password=password)

            if not user:
                form.add_error(None, "Invalid password.")
            elif user and not user.is_active:
                form.add_error(None, "This user may not log in.")
            elif user and user.is_active:
                auth.login(request, user)
                logger.info(f"logged in user.id={user.id}, user.email={user.email}")
                messages.success(request, "Login successful!")
                return redirect(next)
            else:
                # This should not happen normally.
                form.add_error(None, "Invalid form processing.")
    else:

        next = request.GET.get('next', '/')
        initial = dict(next=next)
        form = LoginForm(initial=initial)

    context = dict(form=form, steps=steps, next=next)
    return render(request, "accounts/login.html", context=context)


def password_reset(request):
    steps = breadcrumb_builder([HOME_ICON, LOGIN_ICON])
    context = dict(steps)

    return auth_views.password_reset(request, extra_context=context,
                                     template_name="accounts/password_reset.html",
                                     subject_template_name="accounts/password_reset_subject.txt",
                                     email_template_name="accounts/password_reset_email.html"
                                     )

def password_reset_done(request):
    steps = breadcrumb_builder([HOME_ICON, LOGIN_ICON])
    context = dict(steps=steps)

    return auth_views.password_reset_done(request, extra_context=context,
                                          template_name="accounts/password_reset_done.html")

def password_reset_confirm(request):
    steps = breadcrumb_builder([HOME_ICON, LOGIN_ICON])
    context = dict(steps=steps)

    return auth_views.password_reset_confirm(request, extra_context=context)

def password_reset_complete(request):
    steps = breadcrumb_builder([HOME_ICON, LOGIN_ICON])
    context = dict(steps=steps)

    return auth_views.password_reset_complete(request, extra_context=context)
