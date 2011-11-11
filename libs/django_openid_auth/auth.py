# django-openid-auth -  OpenID integration for django.contrib.auth
#
# Copyright (C) 2008-2010 Canonical Ltd.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Glue between OpenID and django.contrib.auth."""

__metaclass__ = type

from django.conf import settings
from django.contrib.auth.models import User, Group
from openid.consumer.consumer import SUCCESS
from openid.extensions import ax, sreg

from django_openid_auth import teams
from django_openid_auth.models import UserOpenID
from django.contrib import messages
from urlparse import urlparse

class IdentityAlreadyClaimed(Exception):
    pass


class OpenIDBackend:
    """A django.contrib.auth backend that authenticates the user based on
    an OpenID response."""

    def get_user(self, user_id):
        try:
            return User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return None

    def authenticate(self, **kwargs):
        """Authenticate the user based on an OpenID response."""
        # Require that the OpenID response be passed in as a keyword
        # argument, to make sure we don't match the username/password
        # calling conventions of authenticate.

        openid_response = kwargs.get('openid_response')
        if openid_response is None:
            return None

        if openid_response.status != SUCCESS:
            return None

        user = None
        try:
            user_openid = UserOpenID.objects.get(
                claimed_id__exact=openid_response.identity_url)
        except UserOpenID.DoesNotExist:
            if getattr(settings, 'OPENID_CREATE_USERS', False):
                user = self.create_user_from_openid(openid_response)
        else:
            user = user_openid.user

        if user is None:
            return None

        if getattr(settings, 'OPENID_UPDATE_DETAILS_FROM_SREG', False):
            details = self._extract_user_details(openid_response)
            self.update_user_details(user, details)

        teams_response = teams.TeamsResponse.fromSuccessResponse(
            openid_response)
        if teams_response:
            self.update_groups_from_teams(user, teams_response)
            self.update_staff_status_from_teams(user, teams_response)

        return user

    def _extract_user_details(self, openid_response):
        email = fullname = first_name = last_name = nickname = None
        sreg_response = sreg.SRegResponse.fromSuccessResponse(openid_response)
        if sreg_response:
            email = sreg_response.get('email')
            fullname = sreg_response.get('fullname')
            nickname = sreg_response.get('nickname')

        # If any attributes are provided via Attribute Exchange, use
        # them in preference.
        fetch_response = ax.FetchResponse.fromSuccessResponse(openid_response)
        if fetch_response:
            # The myOpenID provider advertises AX support, but uses
            # attribute names from an obsolete draft of the
            # specification.  We check for them first so the common
            # names take precedence.
            email = fetch_response.getSingle(
                'http://schema.openid.net/contact/email', email)
            fullname = fetch_response.getSingle(
                'http://schema.openid.net/namePerson', fullname)
            nickname = fetch_response.getSingle(
                'http://schema.openid.net/namePerson/friendly', nickname)

            email = fetch_response.getSingle(
                'http://axschema.org/contact/email', email)
            fullname = fetch_response.getSingle(
                'http://axschema.org/namePerson', fullname)
            first_name = fetch_response.getSingle(
                'http://axschema.org/namePerson/first', first_name)
            last_name = fetch_response.getSingle(
                'http://axschema.org/namePerson/last', last_name)
            nickname = fetch_response.getSingle(
                'http://axschema.org/namePerson/friendly', nickname)

        if fullname and not (first_name or last_name):
            # Django wants to store first and last names separately,
            # so we do our best to split the full name.
            if ' ' in fullname:
                first_name, last_name = fullname.rsplit(None, 1)
            else:
                first_name = u''
                last_name = fullname

        return dict(email=email, nickname=nickname,
                    first_name=first_name, last_name=last_name)

    def create_user_from_openid(self, openid_response):
        details = self._extract_user_details(openid_response)
        nickname = details['nickname'] or 'openiduser'
        email = details['email'] or ''

        # Pick a username for the user based on their nickname,
        # checking for conflicts.
        i = 1
        while True:
            username = nickname
            if i > 1:
                username += str(i)
            try:
                User.objects.get(username__exact=username)
            except User.DoesNotExist:
                break
            i += 1

        user = User.objects.create_user(username, email, password=None)
        self.update_user_details(user, details)

        user = self.associate_openid(user, openid_response)
        return user

    def associate_openid(self, user, openid_response):
        """Associate an OpenID with a user account."""
        # Check to see if this OpenID has already been claimed.
        try:
            user_openid = UserOpenID.objects.get(
                claimed_id__exact=openid_response.identity_url)
        except UserOpenID.DoesNotExist:
            
            # we need to merge by the email here,
            # only allow a single merge per user
            others = User.objects.filter(email=user.email, profile__openid_merge=False)
            
            # will only trust certain OpenID providers to do the right thing
            # otherwise one can hijack accounts via phony OpenID providers
            url = urlparse(openid_response.identity_url)
            
            # let us know if you need more providers
            trusted = False
            for provider in ( 'www.google.com',  'me.yahoo.com', 'myopenid.com',
                            'livejournal.com', 'blogspot.com', 'openid.aol.com',
                            'wordpress.com'):
                trusted = trusted or url.netloc.endswith(provider)
            
            # you can override migration from the settings
            trusted = trusted and settings.ALLOW_MIGRATION                       
            
            # this merges the authenticated user with an existing user
            if others and trusted:
                print '*** merging an existing user into this openid'
                print '*** openid url %s' % openid_response.identity_url
                # delete the old user
                user.delete()
                
                #  replace with the new user and update the profile accordingly
                user = others[0]
                user.save()
                
                # update the profile with the new information
                user.profile.openid_merge = True
                user.profile.openid = openid_response.identity_url
                user.profile.save()
                

            user_openid = UserOpenID(
                user=user,
                claimed_id=openid_response.identity_url,
                display_id=openid_response.endpoint.getDisplayIdentifier())
            user_openid.save()
        else:
            if user_openid.user != user:
                raise IdentityAlreadyClaimed(
                    "The identity %s has already been claimed"
                    % openid_response.identity_url)

        return user

    def update_user_details(self, user, details):
        
        updated = False
        if details['first_name']:
            user.first_name = details['first_name']
            updated = True

        if details['last_name']:
            user.last_name = details['last_name']
            updated = True

        if details['email']:
            user.email = details['email']
            updated = True

        if not user.get_full_name():
            user.first_name = 'Openid'
            user.last_name  = 'User'
            updated = True

        if updated:
            user.save()

    def update_groups_from_teams(self, user, teams_response):
        teams_mapping_auto = getattr(settings, 'OPENID_LAUNCHPAD_TEAMS_MAPPING_AUTO', False)
        teams_mapping_auto_blacklist = getattr(settings, 'OPENID_LAUNCHPAD_TEAMS_MAPPING_AUTO_BLACKLIST', [])
        teams_mapping = getattr(settings, 'OPENID_LAUNCHPAD_TEAMS_MAPPING', {})
        if teams_mapping_auto:
            #ignore teams_mapping. use all django-groups
            teams_mapping = dict()
            all_groups = Group.objects.exclude(name__in=teams_mapping_auto_blacklist)
            for group in all_groups:
                teams_mapping[group.name] = group.name

        if len(teams_mapping) == 0:
            return

        current_groups = set(user.groups.filter(
                name__in=teams_mapping.values()))
        desired_groups = set(Group.objects.filter(
                name__in=[teams_mapping[lp_team]
                          for lp_team in teams_response.is_member
                          if lp_team in teams_mapping]))
        for group in current_groups - desired_groups:
            user.groups.remove(group)
        for group in desired_groups - current_groups:
            user.groups.add(group)

    def update_staff_status_from_teams(self, user, teams_response):
        if not hasattr(settings, 'OPENID_LAUNCHPAD_STAFF_TEAMS'):
            return

        staff_teams = getattr(settings, 'OPENID_LAUNCHPAD_STAFF_TEAMS', [])
        user.is_staff = False

        for lp_team in teams_response.is_member:
            if lp_team in staff_teams:
                user.is_staff = True
                break

        user.save()

