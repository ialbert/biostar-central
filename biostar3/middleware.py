from __future__ import absolute_import, division, print_function, unicode_literals
import logging

from django.contrib import messages
from django.conf import settings
from django.contrib.auth import get_user_model
from django.http import HttpResponsePermanentRedirect as Redirect
from biostar3.forum.models import UserGroup, GroupPerm
from django.contrib.sites.models import Site
from allauth.socialaccount.adapter import DefaultSocialAccountAdapter
from django.core.cache import cache
from biostar3.forum import cache

# Get current site
User = get_user_model()

logger = logging.getLogger("biostar")

class AutoSignupAdapter(DefaultSocialAccountAdapter):
    def pre_social_login(self, request, sociallogin):

        # This social login already exists.
        if sociallogin.is_existing:
            return

        try:

            # The provider that produces the login
            provider_id = sociallogin.account.get_provider().id

            # Check if we could/should connect it.
            email = sociallogin.account.extra_data.get('email')

            logger.info("connecting %s with %s" % (email, provider_id))

            # Try to get the verification for the account
            verified = sociallogin.account.extra_data.get('verified_email')

            # We will trust some social account providers with the email information.
            verified = verified or (provider_id in settings.TRUSTED_SOCIALACCOUNT_PROVIDERS)

            if email:
                user = User.objects.get(email=email)
                if verified:
                    sociallogin.connect(request, user)
                else:
                    msg = "Attempt to log with email from non verified provider!"
                    logger.error(msg)
                    raise Exception(msg)

        except User.DoesNotExist:
            pass


class GlobalMiddleware(object):
    """Performs tasks that are applied on every request"""

    def process_request(self, request):

        # Set the group based on subdomain on the current request.
        subdomain = settings.GET_SUBDOMAIN(request)

        # Get the group on each request.
        request.group = cache.get_group(subdomain)

        if not request.group:
            # Unable to find the domain redirect.
            # This can cause endless redirects if
            # the DEFAULT_SUBDOMAINS are not set properly.
            site = Site.objects.get_current()
            url = "%s://%s" % (request.scheme, site.domain)
            return Redirect(url)

        # Ensures that requests have all the information needed.
        user = request.user

        # Group level access rights.
        request.user.is_moderator = False

        if not user.is_authenticated():
            # Add minimal user attributes to anonymous users
            user.id = 0
            user.name = "Anonymous"
        else:
            # Check for group permissions.
            perm = GroupPerm.objects.filter(user=user, usergroup=request.group).first()
            user.is_moderator = bool(perm)