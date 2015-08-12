from __future__ import absolute_import, division, print_function, unicode_literals

import json
import logging

# Python modules.
from . import auth

# Django specific modules.
from django.shortcuts import render, redirect
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login
from ratelimit.decorators import ratelimit
from django.contrib.sites.models import Site
from . import query, models
from biostar3.utils.compat import *

logger = logging.getLogger('biostar')

# Get custom user model.
User = get_user_model()


def user_list(request):
    """
    Generates the user listing.
    """
    template_name = "user_list.html"

    q = request.GET.get('q', '')
    if q:
        users = User.objects.filter(name__icontains=q)
    else:
        users = User.objects.all()

    users = users.select_related("profile")

    paginator = query.ExtendedPaginator(request, object_list=users,
                                        time_class=query.TimeLimitValidator,
                                        sort_class=query.UserSortValidator, per_page=25, orphans=False)
    page = paginator.curr_page()

    context = dict(page=page, users=page.object_list, q=q)
    return render(request, template_name, context)


@auth.valid_user
def user_view(request, pk, target=None):
    """
    Generates a single user view.
    The decorator will set the user parameter.
    """
    template_name = "user_view.html"
    # Latest posts
    posts = models.Post.objects.filter(author=target).select_related("root").order_by("-creation_date")[:10]

    target.editable = (request.user == target)

    target.can_be_moderated = auth.can_moderate_user(request=request, target=target, user=request.user)


    top_count = target.post_count(types=models.Post.TOP_LEVEL)
    answer_count = target.post_count(types=[models.Post.ANSWER])
    comment_count = target.post_count(types=[models.Post.COMMENT])

    context = dict(target=target, posts=posts, top_count=top_count,
                   answer_count=answer_count, comment_count=comment_count)

    return render(request, template_name, context)

@login_required
def me_view(request):
    """
    A shortcut that redirects user to their account.
    """
    return redirect(request.user.get_absolute_url())


def badge_view(request, pk):
    """
    Shows awards of one badge.
    """
    template_name = "badge_view.html"
    badge = models.Badge.objects.filter(pk=pk).first()
    if not badge:
        messages.error(request, "Invalid badge selected")
        return redirect("badge_list")


    awards = models.Award.objects.filter(badge=badge).select_related("user", "post", "user__profile").order_by("-date")

    paginator = query.ExtendedPaginator(request,
                                        object_list=awards, per_page=50)

    page = paginator.curr_page()

    context = dict(badge=badge, page=page, awards=page.object_list)

    return render(request, template_name, context)

def badge_list(request):
    """
    Show all the badges
    """
    template_name = "badge_list.html"
    badges = models.Badge.objects.all().order_by("-count")
    context = dict(badges=badges)
    return render(request, template_name, context)

@auth.valid_user
def award_list(request, pk, target=None):
    """
    Show all the awards for a user
    """
    template_name = "award_list.html"
    awards = models.Award.objects.filter(user=target).select_related("badge", "user", "user__profile").order_by("-date")
    paginator = query.ExtendedPaginator(request,
                                        object_list=awards, per_page=50)

    page = paginator.curr_page()

    context = dict(page=page, awards=page.object_list, target=target)
    return render(request, template_name, context)

@login_required
def edit_my_profile(request):
    """
    A shortcut that redirects user to their profile edit.
    """
    return redirect("user_edit", pk=request.user.id)


from allauth.account.views import LoginView


class Login(LoginView):
    """
    This is required to intercept signups from different subdomains.
    Authentication only works for the same domain.
    """

    def dispatch(self, request, *args, **kwargs):
        # This is not used anymore.
        return super(Login, self).dispatch(request, *args, **kwargs)

@ratelimit(key='ip', rate=settings.SIGNUP_RATELIMIT)
def sign_up(request):
    """
    Log In and Sign up is integrated into the same action.

    Performs a series of checks:

    - If the email/password is correct logs the user in.
    - If email exists but password incorrect asks for retry.
    - If no email exists and user provides a password then creates a new user.
    - If Google Recaptcha is active then requires that to pass.

    There are no nested conditionals. As soon as an error is seen we bail out.
    """
    email = request.POST.get('email', '').strip()
    password = request.POST.get('password', '')
    signup = request.POST.get('signup', '')

    recaptcha_response = request.POST.get('g-recaptcha-response')

    login_redirect = redirect(reverse("account_login"))

    was_limited = getattr(request, 'limited', False)

    if was_limited:
        # Rate limits apply to signups.
        messages.error(request, "Access denied! Too many login attempts from the same IP address. Please try later!")
        return login_redirect

    if request.method == 'GET':
        # Get requests go to login page.
        return login_redirect

    if not (email and password):
        # The form requires both email and password.
        messages.error(request, "Form requires an email and a password!")
        return login_redirect

    if not auth.valid_captcha(request):
        # Captcha validation failed.
        return login_redirect

    # At this point the parameters are correct, try to authenticate
    user = authenticate(username=email, password=password)
    if user is not None:
        if user.is_active:
            login(request, user)
            messages.success(request, "Succesfully signed in as %s" % user.name)
            return redirect(reverse("home"))
        else:
            messages.error(request, "This user account has been disabled")
            return login_redirect

    if User.objects.filter(email=email).first():
        # See if the email actually exists.
        messages.error(request, "Incorrect user password.")
        return login_redirect

    if not signup:
        # The user has not requested a sigup.
        messages.error(request, "The email does not exist. If you want to use this email to create an account check the signup checkbox.")
        return login_redirect

    # Try to sign up the user
    user = User.objects.filter(email=email).first()
    if user:
        # The email exists so the password must not match
        messages.error(request, "User password is not correct")
        return login_redirect

    try:
        # Now sign up and log in the user.
        user = User.objects.create(email=email)
        user.set_password(password)
        user.save()
        user = authenticate(username=email, password=password)
        login(request, user)
        messages.success(request, "Succesfully signed up as %s" % user.name)
        return redirect(reverse("home"))
    except Exception as exc:
        logger.error(exc)
        messages.error(request, "Unable to create user")
        return login_redirect

    messages.error(request, "Site login error. Please contact the admins.")
    return redirect(reverse("account_login"))


