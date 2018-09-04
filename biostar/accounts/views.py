import uuid, logging

from django.contrib import messages
from ratelimit.decorators import ratelimit
from django.shortcuts import render, redirect
from django.urls import reverse
from django.contrib.auth import views as auth_views
from django.utils.safestring import mark_safe
from django.conf import settings
from django.contrib.auth import logout, login
from django.contrib.auth.decorators import login_required
from django.utils.encoding import force_text, force_bytes
from django.utils.http import urlsafe_base64_decode

from .tokens import account_verification_token
from . import forms
from .models import User, Profile
from .auth import check_user, send_verification_email
from .util import now
from .const import *

logger = logging.getLogger('engine')


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def edit_profile(request):

    if request.user.is_anonymous:
        messages.error(request, "Must be logged in to edit profile")
        return redirect("/")

    user = request.user
    form = forms.EditProfile(user=user)

    if request.method == "POST":
        form = forms.EditProfile(data=request.POST, user=user)
        if form.is_valid():
            form.save(request=request)
            return redirect(reverse("public_profile", kwargs=dict(uid=user.profile.uid)))

    context = dict(user=user, form=form)
    return render(request, 'accounts/edit_profile.html', context)


@login_required
def toggle_moderate(request):

    user = request.user

    if settings.ALLOW_SELF_MODERATE:

        role = Profile.NORMAL if user.profile.is_moderator else Profile.MODERATOR
        Profile.objects.filter(user=user).update(role=role)
        mapper = {Profile.MODERATOR:" a moderator"}
        messages.success(request, f"You are now {mapper.get(role, 'not a moderator')}")

    return redirect(reverse("public_profile", kwargs=dict(uid=user.profile.uid)))


@login_required
def user_moderate(request, uid):

    source = request.user
    target = User.objects.filter(profile__uid=uid).first()
    form = forms.UserModerate(source=source, target=target, request=request)

    if request.method == "POST":

        form = forms.UserModerate(source=source, data=request.POST, target=target, request=request)
        if form.is_valid():
            form.save()
            return redirect(reverse("public_profile", kwargs=dict(uid=uid)))
        else:
            msg = ','.join([y for x in form.errors.values() for y in x])
            messages.error(request, msg)

    context = dict(form=form, target=target)
    return render(request, "accounts/user_moderate.html", context)


def public_profile(request, uid):

    user_profile = Profile.objects.filter(uid=uid).first()
    forum_enabled = settings.ENABLE_FORUM or settings.ONLY_ENABLE_FORUM

    # Get the active tab, defaults to project
    active_tab = request.GET.get(ACTIVE_TAB, HAS_POSTS)

    if not user_profile:
        messages.error(request, "User profile does not exist")
        return redirect("/")

    active_tab = active_tab if (active_tab in PROFILE_TABS) else HAS_POSTS
    can_moderate = user_profile.can_moderate(source=request.user)

    context = dict(user=user_profile.user, enable_forum=forum_enabled,
                   const_name=ACTIVE_TAB, const_post=HAS_POSTS, const_project=HAS_PROJECT,
                   const_recipes=HAS_RECIPES, can_moderate=can_moderate)

    context.update({active_tab: ACTIVE_TAB})

    return render(request, 'accounts/public_profile.html', context)


@login_required
def profile(request):
    uid = request.user.profile.uid

    return public_profile(request, uid)


def toggle_notify(request):

    if request.user.is_anonymous:
        messages.error(request, "Must be logged in to edit profile")
        return redirect("/")

    user = request.user
    user.profile.notify = not user.profile.notify
    user.profile.save()

    msg = "Emails notifications disabled."
    if user.profile.notify:
        msg = "Emails notifications enabled."

    messages.success(request, msg)
    return redirect(reverse('public_profile', kwargs=dict(uid=user.profile.uid)))


@ratelimit(key='ip', rate='10/m', block=True, method=ratelimit.UNSAFE)
def user_signup(request):

    if not settings.ALLOW_SIGNUP:
        messages.info(request, "Signups are not enabled on this site.")
        return redirect("/")

    if request.method == 'POST':

        form = forms.SignUpWithCaptcha(request.POST)
        if form.is_valid():
            email = form.cleaned_data.get('email')
            password = form.cleaned_data.get('password1')
            name = email.split("@")[0]
            user = User.objects.create(email=email, first_name=name)
            user.set_password(password)
            user.username = name.split()[0] + str(user.id)
            user.save()

            #login(request, user)
            #Profile.objects.filter(user=user).update(last_login=now())
            #messages.success(request, "Login successful!")

            send_verification_email(user=user)
            logger.info(f"Signed up user.id={user.id}, user.email={user.email}")
            msg = mark_safe("Signup successful! <b>Please verify your email to complete registration.</b>")
            messages.info(request, msg)

            return redirect("/")
    else:
        form = forms.SignUpWithCaptcha()
    context = dict(form=form, captcha_site_key=settings.RECAPTCHA_PUBLIC_KEY)
    return render(request, 'accounts/signup.html', context=context)


def user_logout(request):

    if request.method == "POST":

        form = forms.LogoutForm(request.POST)

        if form.is_valid():
            logout(request)
            messages.info(request, "You have been logged out")
            return redirect("/")

    form = forms.LogoutForm()

    context = dict(form=form)

    return render(request, "accounts/logout.html", context=context)


def user_login(request):

    form = forms.LoginForm()
    if request.method == "POST":
        form = forms.LoginForm(data=request.POST)

        if form.is_valid():

            email = form.cleaned_data['email']
            password = form.cleaned_data['password']

            user = User.objects.filter(email__iexact=email).order_by('-id').first()
            message, valid_user = check_user(email=email, password=password)

            if valid_user:
                login(request, user)
                Profile.objects.filter(user=user).update(last_login=now())
                messages.success(request, "Login successful!")
                return redirect("/")
            else:
                messages.error(request, mark_safe(message))

        messages.error(request, mark_safe(form.errors))

    context = dict(form=form)
    return render(request, "accounts/login.html", context=context)


@login_required
def send_email_verify(request):
    "Send one-time valid link to validate email"

    # Sends verification email with a token
    user = request.user

    send_verification_email(user=user)

    messages.success(request, "Verification sent, check your email.")

    return redirect(reverse("public_profile", kwargs=dict(uid=user.profile.uid)))


def email_verify_account(request, uidb64, token):
    "Verify one time link sent to a users email"

    uid = force_text(urlsafe_base64_decode(uidb64))
    user = User.objects.filter(pk=uid).first()

    if user and account_verification_token.check_token(user, token):
        Profile.objects.filter(user=user).update(email_verified=True)
        login(request, user)
        messages.success(request, "Email verified!")
        return redirect(reverse('public_profile', kwargs=dict(uid=user.profile.uid)))

    messages.error(request, "Link is expired.")
    return redirect("/")


def password_reset(request):
    context = dict()

    return auth_views.password_reset(request, extra_context=context,
                                     template_name="accounts/password_reset_form.html",
                                     subject_template_name="accounts/password_reset_subject.txt",
                                     email_template_name="accounts/password_reset_email.html"
                                     )

def password_reset_done(request):
    context = dict()
    return auth_views.password_reset_done(request, extra_context=context,
                                          template_name="accounts/password_reset_done.html")


def pass_reset_confirm(request, uidb64, token):
    context = dict()

    return auth_views.password_reset_confirm(request, extra_context=context,
                                             template_name="accounts/password_reset_confirm.html",
                                            uidb64=uidb64, token=token)

def password_reset_complete(request):
    context = dict()

    return auth_views.password_reset_complete(request, extra_context=context,
                                              template_name="accounts/password_reset_complete.html",)
