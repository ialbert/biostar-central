from __future__ import absolute_import, division, print_function, unicode_literals

import urllib, urllib2, json, logging

# Python modules.
from collections import OrderedDict, defaultdict
from . import auth

# Django specific modules.
from django.shortcuts import render, redirect
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login
from ratelimit.decorators import ratelimit, is_ratelimited
from django.contrib.sites.models import Site
from . import query, models

logger = logging.getLogger('biostar')

# Get custom user model.
User = get_user_model()

# Loads this group when none are specified.
DEFAULT_GROUP = models.UserGroup.objects.filter(name=settings.DEFAULT_GROUP_NAME).first()

def user_list(request):
    """
    Generates the user listing.
    """
    template_name = "user_list.html"

    q = request.GET.get('q','')
    if q:
        users = User.objects.filter(name__icontains=q)
    else:
        users = User.objects.all()

    paginator = query.ExtendedPaginator(request, object_list=users,
                                         time_class=query.TimeLimitValidator,
                                        sort_class=query.UserSortValidator, per_page=25, orphans=False)
    page = paginator.curr_page()

    context = dict(page=page, users=page.object_list, q=q)
    return render(request, template_name, context)

@auth.valid_user
def user_view(request, pk, user=None):
    """
    Generates a single user view.
    The decorator will set the user parameter.
    """

    template_name = "user_view.html"
    context = dict(target=user)
    return render(request, template_name, context)

@login_required
def me_view(request):
    """
    A shortcut that redirects user to their account.
    """
    return request.user.get_absolute_url()


from allauth.account.views import LoginView
class Login(LoginView):
    """
    This is required to intercept signups from different subdomains.
    Authentication only works for the same domain.
    """
    def dispatch(self, request, *args, **kwargs):
        # Override the next parameter if the user is visiting a different group.
        # This is necessary as OAuth will only work for a single domain.
        # After log in the user is redirected to the group site.
        group = request.group
        if group.domain != DEFAULT_GROUP.domain:
            site = Site.objects.get_current()
            login_url = reverse("account_login")
            next_url = reverse("group_redirect", kwargs=dict(domain=group.domain))
            params = dict(
                scheme=request.scheme,
                next_url=next_url,
                domain=site.domain,
                login_url=login_url,
            )
            site_url = "%(scheme)s://%(domain)s%(login_url)s?next=%(next_url)s" % params
            return redirect(site_url)

        return super(Login, self).dispatch(request, *args, **kwargs)

CAPTCHA_VERIFY_URL = "https://www.google.com/recaptcha/api/siteverify"
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

    if settings.RECAPTCHA_PUBLIC_KEY and not recaptcha_response:
        # The capthca is active but the user has not solved it.
        messages.error(request, "Please solve the captcha!")
        return login_redirect

    if settings.RECAPTCHA_PUBLIC_KEY:
        # Verfiying the captcha response.

        data = {
            'secret': settings.RECAPTCHA_SECRET_KEY,
            'response': recaptcha_response,
            'remoteip': auth.remote_ip(request),
        }

        try:
            # Validate the captcha.
            data = urllib.urlencode(data)
            conn = urllib2.Request(CAPTCHA_VERIFY_URL, data)
            response = urllib2.urlopen(conn).read()
            result = json.loads(response)
            if not result.get('success'):
                # User failed at solving the capthca.
                messages.error(request, "Failed at the captcha authentication. Please try again!")
                return login_redirect
        except Exception, exc:
            # This here is triggered on unexpected errors while solving the capthca.
            logger.error(exc)
            messages.error(request, "Unable to complete captcha challenge: %s" % exc)
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
        messages.error(request, "If you want to create an account check the signup checkbox.")
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
    except Exception, exc:
        logger.error(exc)
        messages.error(request, "Unable to create user")
        return login_redirect

    messages.error(request, "Site login error. Please contact the admins.")
    return redirect(reverse("account_login"))


