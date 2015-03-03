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

logger = logging.getLogger('biostar')

# Get custom user model.
User = get_user_model()

@login_required
def me_view(request):
    """
    A shortcut that redirects user to their account.
    """
    return request.user.get_absolute_url()


CAPTCHA_VERIFY_URL = "https://www.google.com/recaptcha/api/siteverify"

def sign_up(request):

    email = request.POST.get('email', '').strip()
    password = request.POST.get('password')
    recaptcha_response = request.POST.get('g-recaptcha-response')

    login_redirect = redirect(reverse("account_login"))

    if request.method == 'GET':
        return login_redirect

    if not email and password:
        messages.error(request, "Please enter an email and a password!")
        return login_redirect

    if settings.RECAPTCHA_PUBLIC_KEY and not recaptcha_response:
        messages.error(request, "Please solve the captcha!")
        return login_redirect

    if settings.RECAPTCHA_PUBLIC_KEY:
        # Verfiying the captcha response.

        data = {
            'secret' : settings.RECAPTCHA_PRIVATE_KEY,
            'response' : recaptcha_response,
            'remoteip' : auth.remote_ip(request),
        }

        try:
            data = urllib.urlencode(data)
            conn = urllib2.Request(CAPTCHA_VERIFY_URL, data)
            response = urllib2.urlopen(conn).read()
            result = json.loads(response)
            if not result.get('success'):
                messages.error(request, "Failed at the captcha authentication. Please try again!")
                return login_redirect
        except Exception, exc:
            logger.error(exc)
            messages.error(request, "Unable to complete captcha challenge: %s" % exc)
            return login_redirect


    # At this point the parameters are correct, try to authenticate
    user = authenticate(username=email, password=password)
    if user is not None:
        if user.is_active:
            login(request, user)
            return redirect(reverse("home"))
        else:
            messages.error(request, "This user account has been disabled")
            return login_redirect

    # If we are here then authentication was not successful.
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
        return redirect(reverse("home"))
    except Exception, exc:
        logger.error(exc)
        messages.error(request, "Unable to create user")
        return login_redirect

    messages.error(request, "Site login error. Please contact the admins.")
    return redirect(reverse("account_login"))


