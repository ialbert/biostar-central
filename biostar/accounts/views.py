import uuid, logging

from django.contrib import messages
from ratelimit.decorators import ratelimit
from django.shortcuts import render, redirect
from django.contrib.auth import views as auth_views
from django.utils.safestring import mark_safe
from django.conf import settings
from django.contrib.auth import logout, login
from django.contrib.auth.decorators import login_required
from django.utils.encoding import force_text, force_bytes
from django.utils.http import urlsafe_base64_decode
from django.contrib.auth.decorators import user_passes_test
from django.core import signing
from django.core.paginator import Paginator

from biostar.utils.decorators import object_exists
from biostar.utils.shortcuts import reverse
from allauth.socialaccount.models import SocialApp
from .tokens import account_verification_token
from . import forms
from .models import User, Profile
from .auth import check_user, send_verification_email, query_topic
from .util import now, get_uuid
from .const import *

logger = logging.getLogger('engine')


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

    profile = Profile.objects.filter(uid=uid).first()

    if not profile:
        messages.error(request, "User does not exist")
        return redirect("/")

    # Get the active tab, defaults to project
    active_tab = request.GET.get(ACTIVE_TAB, POSTS)

    can_moderate = profile.can_moderate(source=request.user)

    objs = query_topic(user=profile.user, request=request, active=active_tab)

    page = request.GET.get("page", 1)

    paginator = Paginator(objs, per_page=20)
    page = page if page is not None else 1

    objs = paginator.get_page(page)

    context = dict(user=profile.user, objs=objs,
                   const_name=ACTIVE_TAB, const_post=POSTS, const_project=PROJECT,
                   const_recipes=RECIPES, can_moderate=can_moderate)

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
            user = form.save()
            login(request, user,  backend="django.contrib.auth.backends.ModelBackend")
            Profile.objects.filter(user=user).update(last_login=now())
            messages.success(request, "Login successful!")
            msg = mark_safe("Signup successful!")
            messages.info(request, msg)

            return redirect("/")
    else:
        form = forms.SignUpWithCaptcha()

    context = dict(form=form, captcha_site_key=settings.RECAPTCHA_PUBLIC_KEY,
                   social_login=SocialApp.objects.all())
    return render(request, 'accounts/signup.html', context=context)


@user_passes_test(lambda u: u.is_superuser)
def debug_user(request):
    """
    Allows superusers to log in as a regular user to troubleshoot problems.
    """

    target = request.GET.get("uid", "")
    profile = Profile.objects.filter(uid=target).first()

    if not profile:
        messages.error(request, "User does not exists.")
        return redirect("/")

    user = profile.user
    login(request, user, backend="django.contrib.auth.backends.ModelBackend")
    messages.success(request, "Login successful!")

    logger.info(f"""uid={request.user.profile.uid} impersonated 
                    uid={profile.uid}.""")

    return redirect("/")


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
                login(request, user,  backend="django.contrib.auth.backends.ModelBackend")
                Profile.objects.filter(user=user).update(last_login=now())
                messages.success(request, "Login successful!")
                return redirect("/")
            else:
                messages.error(request, mark_safe(message))

        messages.error(request, mark_safe(form.errors))

    context = dict(form=form, social_login=SocialApp.objects.all())
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
        login(request, user, backend="django.contrib.auth.backends.ModelBackend")
        messages.success(request, "Email verified!")
        return redirect(reverse('public_profile', kwargs=dict(uid=user.profile.uid)))

    messages.error(request, "Link is expired.")
    return redirect("/")


def external_login(request):

    payload = request.GET.get("payload", "")

    try:
        signer = signing.Signer(settings.LOGIN_PRIVATE_KEY)
        user_email = signer.unsign(payload)
        user = User.objects.filter(email=user_email).first()

        if not user:
            name = user_email.split("@")[0]
            user = User.objects.create(email=user_email, first_name=name,
                                       password=str(get_uuid(16)))
            user.username = name.split()[0] + str(get_uuid(8))
            user.save()
            msg = f"Signed up, <a href={reverse('password_reset')}><b> Please reset your password.</b></a>"
            messages.success(request, mark_safe(msg))

        login(request, user, backend="django.contrib.auth.backends.ModelBackend")
        messages.success(request, "Logged in!")
        return redirect("/")

    except Exception as exc:
        logger.error(f"Error:{exc}")
        messages.error(request, "Error unsigning.")

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
