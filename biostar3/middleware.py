from __future__ import absolute_import, division, print_function, unicode_literals

from django.contrib import messages
from django.conf import settings
from django.contrib.auth import get_user_model
from django.http import HttpResponsePermanentRedirect as Redirect
from biostar3.forum.models import Group
from django.contrib.sites.models import Site
from allauth.socialaccount.adapter import DefaultSocialAccountAdapter

# Get current site
User = get_user_model()

SITE = Site.objects.get_current()

# Loads this group when none are specified.
DEFAULT_GROUP = Group.objects.filter(name=settings.DEFAULT_GROUP_NAME).first()

def full_url(request, url):
    if request.is_secure():
        return "https://%s" % url
    else:
        return "http://%s" % url

class AutoSignupAdapter(DefaultSocialAccountAdapter):

    def pre_social_login(self, request, sociallogin):

        # This social login already exists.
        if sociallogin.is_existing:
            return

        try:
            print (sociallogin.account.extra_data)

            # Check if we could/should connect it.
            email = sociallogin.account.extra_data.get('email')
            #verified = sociallogin.account.extra_data.get('verified_email')
            if email:
                user = User.objects.get(email=email)
                sociallogin.connect(request, user)
        except User.DoesNotExist:
            pass

class GlobalMiddleware(object):
    """Performs tasks that are applied on every request"""

    def process_request(self, request):

        # Ensures that requests have all the information needed.
        user = request.user
        if not user.is_authenticated():
            user.is_moderator = user.is_admin = False

        # Set the group based on subdomain on the current request
        subdomain = settings.GET_SUBDOMAIN(request)
        if subdomain in settings.DEFAULT_SUBDOMAINS:
            group = DEFAULT_GROUP
        else:
            group = Group.objects.filter(name__iexact=subdomain).first()
            if not group:
                return Redirect(full_url(request, SITE.domain))

        # Groups need to be set on each request.
        request.group = group

