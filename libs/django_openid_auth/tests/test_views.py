# django-openid-auth -  OpenID integration for django.contrib.auth
#
# Copyright (C) 2009-2010 Canonical Ltd.
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

import cgi
import unittest

from django.conf import settings
from django.contrib.auth.models import User, Group
from django.http import HttpRequest
from django.test import TestCase
from openid.consumer.consumer import SuccessResponse
from openid.extensions import ax, sreg
from openid.fetchers import (
    HTTPFetcher, HTTPFetchingError, HTTPResponse, setDefaultFetcher)
from openid.oidutil import importElementTree
from openid.server.server import BROWSER_REQUEST_MODES, ENCODE_URL, Server
from openid.store.memstore import MemoryStore

from django_openid_auth import teams
from django_openid_auth.models import UserOpenID
from django_openid_auth.views import sanitise_redirect_url
from django_openid_auth.signals import openid_login_complete


ET = importElementTree()


class StubOpenIDProvider(HTTPFetcher):

    def __init__(self, base_url):
        self.store = MemoryStore()
        self.identity_url = base_url + 'identity'
        self.localid_url = base_url + 'localid'
        self.endpoint_url = base_url + 'endpoint'
        self.server = Server(self.store, self.endpoint_url)
        self.last_request = None
        self.type_uris = ['http://specs.openid.net/auth/2.0/signon']

    def fetch(self, url, body=None, headers=None):
        if url == self.identity_url:
            # Serve an XRDS document directly, pointing at our endpoint.
            type_uris = ['<Type>%s</Type>' % uri for uri in self.type_uris]
            return HTTPResponse(
                url, 200, {'content-type': 'application/xrds+xml'}, """\
<?xml version="1.0"?>
<xrds:XRDS
    xmlns="xri://$xrd*($v*2.0)"
    xmlns:xrds="xri://$xrds">
  <XRD>
    <Service priority="0">
      %s
      <URI>%s</URI>
      <LocalID>%s</LocalID>
    </Service>
  </XRD>
</xrds:XRDS>
""" % ('\n'.join(type_uris), self.endpoint_url, self.localid_url))
        elif url.startswith(self.endpoint_url):
            # Gather query parameters
            query = {}
            if '?' in url:
                query.update(cgi.parse_qsl(url.split('?', 1)[1]))
            if body is not None:
                query.update(cgi.parse_qsl(body))
            self.last_request = self.server.decodeRequest(query)

            # The browser based requests should not be handled through
            # the fetcher interface.
            assert self.last_request.mode not in BROWSER_REQUEST_MODES

            response = self.server.handleRequest(self.last_request)
            webresponse = self.server.encodeResponse(response)
            return HTTPResponse(url,  webresponse.code, webresponse.headers,
                                webresponse.body)
        else:
            raise HTTPFetchingError('unknown URL %s' % url)

    def parseFormPost(self, content):
        """Parse an HTML form post to create an OpenID request."""
        # Hack to make the javascript XML compliant ...
        content = content.replace('i < elements.length',
                                  'i &lt; elements.length')
        tree = ET.XML(content)
        form = tree.find('.//form')
        assert form is not None, 'No form in document'
        assert form.get('action') == self.endpoint_url, (
            'Form posts to %s instead of %s' % (form.get('action'),
                                                self.endpoint_url))
        query = {}
        for input in form.findall('input'):
            if input.get('type') != 'hidden':
                continue
            query[input.get('name').encode('UTF-8')] = \
                input.get('value').encode('UTF-8')
        self.last_request = self.server.decodeRequest(query)
        return self.last_request


class RelyingPartyTests(TestCase):
    urls = 'django_openid_auth.tests.urls'

    def setUp(self):
        super(RelyingPartyTests, self).setUp()
        self.provider = StubOpenIDProvider('http://example.com/')
        setDefaultFetcher(self.provider, wrap_exceptions=False)

        self.old_login_redirect_url = getattr(settings, 'LOGIN_REDIRECT_URL', '/accounts/profile/')
        self.old_create_users = getattr(settings, 'OPENID_CREATE_USERS', False)
        self.old_update_details = getattr(settings, 'OPENID_UPDATE_DETAILS_FROM_SREG', False)
        self.old_sso_server_url = getattr(settings, 'OPENID_SSO_SERVER_URL', None)
        self.old_teams_map = getattr(settings, 'OPENID_LAUNCHPAD_TEAMS_MAPPING', {})
        self.old_use_as_admin_login = getattr(settings, 'OPENID_USE_AS_ADMIN_LOGIN', False)

        settings.OPENID_CREATE_USERS = False
        settings.OPENID_UPDATE_DETAILS_FROM_SREG = False
        settings.OPENID_SSO_SERVER_URL = None
        settings.OPENID_LAUNCHPAD_TEAMS_MAPPING = {}
        settings.OPENID_USE_AS_ADMIN_LOGIN = False

    def tearDown(self):
        settings.LOGIN_REDIRECT_URL = self.old_login_redirect_url
        settings.OPENID_CREATE_USERS = self.old_create_users
        settings.OPENID_UPDATE_DETAILS_FROM_SREG = self.old_update_details
        settings.OPENID_SSO_SERVER_URL = self.old_sso_server_url
        settings.OPENID_LAUNCHPAD_TEAMS_MAPPING = self.old_teams_map
        settings.OPENID_USE_AS_ADMIN_LOGIN = self.old_use_as_admin_login

        setDefaultFetcher(None)
        super(RelyingPartyTests, self).tearDown()

    def complete(self, openid_response):
        """Complete an OpenID authentication request."""
        # The server can generate either a redirect or a form post
        # here.  For simplicity, force generation of a redirect.
        openid_response.whichEncoding = lambda: ENCODE_URL
        webresponse = self.provider.server.encodeResponse(openid_response)
        self.assertEquals(webresponse.code, 302)
        redirect_to = webresponse.headers['location']
        self.assertTrue(redirect_to.startswith(
                'http://testserver/openid/complete/'))
        return self.client.get('/openid/complete/',
            dict(cgi.parse_qsl(redirect_to.split('?', 1)[1])))

    def test_login(self):
        user = User.objects.create_user('someuser', 'someone@example.com')
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()

        # The login form is displayed:
        response = self.client.get('/openid/login/')
        self.assertTemplateUsed(response, 'openid/login.html')

        # Posting in an identity URL begins the authentication request:
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity',
             'next': '/getuser/'})
        self.assertContains(response, 'OpenID transaction in progress')

        openid_request = self.provider.parseFormPost(response.content)
        self.assertEquals(openid_request.mode, 'checkid_setup')
        self.assertTrue(openid_request.return_to.startswith(
                'http://testserver/openid/complete/'))

        # Complete the request.  The user is redirected to the next URL.
        openid_response = openid_request.answer(True)
        response = self.complete(openid_response)
        self.assertRedirects(response, 'http://testserver/getuser/')

        # And they are now logged in:
        response = self.client.get('/getuser/')
        self.assertEquals(response.content, 'someuser')

    def test_login_no_next(self):
        """Logins with no next parameter redirect to LOGIN_REDIRECT_URL."""
        user = User.objects.create_user('someuser', 'someone@example.com')
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()

        settings.LOGIN_REDIRECT_URL = '/getuser/'
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity'})
        self.assertContains(response, 'OpenID transaction in progress')

        openid_request = self.provider.parseFormPost(response.content)
        self.assertEquals(openid_request.mode, 'checkid_setup')
        self.assertTrue(openid_request.return_to.startswith(
                'http://testserver/openid/complete/'))

        # Complete the request.  The user is redirected to the next URL.
        openid_response = openid_request.answer(True)
        response = self.complete(openid_response)
        self.assertRedirects(
            response, 'http://testserver' + settings.LOGIN_REDIRECT_URL)

    def test_login_sso(self):
        settings.OPENID_SSO_SERVER_URL = 'http://example.com/identity'
        user = User.objects.create_user('someuser', 'someone@example.com')
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()

        # Requesting the login form immediately begins an
        # authentication request.
        response = self.client.get('/openid/login/', {'next': '/getuser/'})
        self.assertEquals(response.status_code, 200)
        self.assertContains(response, 'OpenID transaction in progress')

        openid_request = self.provider.parseFormPost(response.content)
        self.assertEquals(openid_request.mode, 'checkid_setup')
        self.assertTrue(openid_request.return_to.startswith(
                'http://testserver/openid/complete/'))

        # Complete the request.  The user is redirected to the next URL.
        openid_response = openid_request.answer(True)
        response = self.complete(openid_response)
        self.assertRedirects(response, 'http://testserver/getuser/')

        # And they are now logged in:
        response = self.client.get('/getuser/')
        self.assertEquals(response.content, 'someuser')

    def test_login_create_users(self):
        settings.OPENID_CREATE_USERS = True
        # Create a user with the same name as we'll pass back via sreg.
        User.objects.create_user('someuser', 'someone@example.com')

        # Posting in an identity URL begins the authentication request:
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity',
             'next': '/getuser/'})
        self.assertContains(response, 'OpenID transaction in progress')

        # Complete the request, passing back some simple registration
        # data.  The user is redirected to the next URL.
        openid_request = self.provider.parseFormPost(response.content)
        sreg_request = sreg.SRegRequest.fromOpenIDRequest(openid_request)
        openid_response = openid_request.answer(True)
        sreg_response = sreg.SRegResponse.extractResponse(
            sreg_request, {'nickname': 'someuser', 'fullname': 'Some User',
                           'email': 'foo@example.com'})
        openid_response.addExtension(sreg_response)
        response = self.complete(openid_response)
        self.assertRedirects(response, 'http://testserver/getuser/')

        # And they are now logged in as a new user (they haven't taken
        # over the existing "someuser" user).
        response = self.client.get('/getuser/')
        self.assertEquals(response.content, 'someuser2')

        # Check the details of the new user.
        user = User.objects.get(username='someuser2')
        self.assertEquals(user.first_name, 'Some')
        self.assertEquals(user.last_name, 'User')
        self.assertEquals(user.email, 'foo@example.com')

    def test_login_update_details(self):
        settings.OPENID_UPDATE_DETAILS_FROM_SREG = True
        user = User.objects.create_user('testuser', 'someone@example.com')
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()

        # Posting in an identity URL begins the authentication request:
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity',
             'next': '/getuser/'})
        self.assertContains(response, 'OpenID transaction in progress')

        # Complete the request, passing back some simple registration
        # data.  The user is redirected to the next URL.
        openid_request = self.provider.parseFormPost(response.content)
        sreg_request = sreg.SRegRequest.fromOpenIDRequest(openid_request)
        openid_response = openid_request.answer(True)
        sreg_response = sreg.SRegResponse.extractResponse(
            sreg_request, {'nickname': 'someuser', 'fullname': 'Some User',
                           'email': 'foo@example.com'})
        openid_response.addExtension(sreg_response)
        response = self.complete(openid_response)
        self.assertRedirects(response, 'http://testserver/getuser/')

        # And they are now logged in as testuser (the passed in
        # nickname has not caused the username to change).
        response = self.client.get('/getuser/')
        self.assertEquals(response.content, 'testuser')

        # The user's full name and email have been updated.
        user = User.objects.get(username='testuser')
        self.assertEquals(user.first_name, 'Some')
        self.assertEquals(user.last_name, 'User')
        self.assertEquals(user.email, 'foo@example.com')

    def test_login_uses_sreg_extra_fields(self):
        # The configurable sreg attributes are used in the request.
        settings.OPENID_SREG_EXTRA_FIELDS = ('language',)
        user = User.objects.create_user('testuser', 'someone@example.com')
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()

        # Posting in an identity URL begins the authentication request:
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity',
             'next': '/getuser/'})

        openid_request = self.provider.parseFormPost(response.content)
        sreg_request = sreg.SRegRequest.fromOpenIDRequest(openid_request)
        for field in ('email', 'fullname', 'nickname', 'language'):
            self.assertTrue(field in sreg_request)

    def test_login_attribute_exchange(self):
        settings.OPENID_UPDATE_DETAILS_FROM_SREG = True
        user = User.objects.create_user('testuser', 'someone@example.com')
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()

        # Configure the provider to advertise attribute exchange
        # protocol and start the authentication process:
        self.provider.type_uris.append('http://openid.net/srv/ax/1.0')
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity',
             'next': '/getuser/'})
        self.assertContains(response, 'OpenID transaction in progress')

        # The resulting OpenID request uses the Attribute Exchange
        # extension rather than the Simple Registration extension.
        openid_request = self.provider.parseFormPost(response.content)
        sreg_request = sreg.SRegRequest.fromOpenIDRequest(openid_request)
        self.assertEqual(sreg_request.required, [])
        self.assertEqual(sreg_request.optional, [])

        fetch_request = ax.FetchRequest.fromOpenIDRequest(openid_request)
        self.assertTrue(fetch_request.has_key(
                'http://axschema.org/contact/email'))
        self.assertTrue(fetch_request.has_key(
                'http://axschema.org/namePerson'))
        self.assertTrue(fetch_request.has_key(
                'http://axschema.org/namePerson/first'))
        self.assertTrue(fetch_request.has_key(
                'http://axschema.org/namePerson/last'))
        self.assertTrue(fetch_request.has_key(
                'http://axschema.org/namePerson/friendly'))
        # myOpenID compatibilty attributes:
        self.assertTrue(fetch_request.has_key(
                'http://schema.openid.net/contact/email'))
        self.assertTrue(fetch_request.has_key(
                'http://schema.openid.net/namePerson'))
        self.assertTrue(fetch_request.has_key(
                'http://schema.openid.net/namePerson/friendly'))

        # Build up a response including AX data.
        openid_response = openid_request.answer(True)
        fetch_response = ax.FetchResponse(fetch_request)
        fetch_response.addValue(
            'http://axschema.org/contact/email', 'foo@example.com')
        fetch_response.addValue(
            'http://axschema.org/namePerson/first', 'Firstname')
        fetch_response.addValue(
            'http://axschema.org/namePerson/last', 'Lastname')
        fetch_response.addValue(
            'http://axschema.org/namePerson/friendly', 'someuser')
        openid_response.addExtension(fetch_response)
        response = self.complete(openid_response)
        self.assertRedirects(response, 'http://testserver/getuser/')

        # And they are now logged in as testuser (the passed in
        # nickname has not caused the username to change).
        response = self.client.get('/getuser/')
        self.assertEquals(response.content, 'testuser')

        # The user's full name and email have been updated.
        user = User.objects.get(username='testuser')
        self.assertEquals(user.first_name, 'Firstname')
        self.assertEquals(user.last_name, 'Lastname')
        self.assertEquals(user.email, 'foo@example.com')

    def test_login_teams(self):
        settings.OPENID_LAUNCHPAD_TEAMS_MAPPING = {'teamname': 'groupname',
                                                   'otherteam': 'othergroup'}
        user = User.objects.create_user('testuser', 'someone@example.com')
        group = Group(name='groupname')
        group.save()
        ogroup = Group(name='othergroup')
        ogroup.save()
        user.groups.add(ogroup)
        user.save()
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()

        # Posting in an identity URL begins the authentication request:
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity',
             'next': '/getuser/'})
        self.assertContains(response, 'OpenID transaction in progress')

        # Complete the request
        openid_request = self.provider.parseFormPost(response.content)
        openid_response = openid_request.answer(True)
        teams_request = teams.TeamsRequest.fromOpenIDRequest(openid_request)
        teams_response = teams.TeamsResponse.extractResponse(
            teams_request, 'teamname,some-other-team')
        openid_response.addExtension(teams_response)
        response = self.complete(openid_response)
        self.assertRedirects(response, 'http://testserver/getuser/')

        # And they are now logged in as testuser
        response = self.client.get('/getuser/')
        self.assertEquals(response.content, 'testuser')

        # The user's groups have been updated.
        user = User.objects.get(username='testuser')
        self.assertTrue(group in user.groups.all())
        self.assertTrue(ogroup not in user.groups.all())

    def test_login_teams_automapping(self):
        settings.OPENID_LAUNCHPAD_TEAMS_MAPPING = {'teamname': 'groupname',
                                                   'otherteam': 'othergroup'}
        settings.OPENID_LAUNCHPAD_TEAMS_MAPPING_AUTO = True
        settings.OPENID_LAUNCHPAD_TEAMS_MAPPING_AUTO_BLACKLIST = ['django-group1', 'django-group2']
        user = User.objects.create_user('testuser', 'someone@example.com')
        group1 = Group(name='django-group1')
        group1.save()
        group2 = Group(name='django-group2')
        group2.save()
        group3 = Group(name='django-group3')
        group3.save()
        user.save()
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()

        # Posting in an identity URL begins the authentication request:
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity',
             'next': '/getuser/'})
        self.assertContains(response, 'OpenID transaction in progress')

        # Complete the request
        openid_request = self.provider.parseFormPost(response.content)
        openid_response = openid_request.answer(True)
        teams_request = teams.TeamsRequest.fromOpenIDRequest(openid_request)

        self.assertEqual(group1 in user.groups.all(), False)
        self.assertEqual(group2 in user.groups.all(), False)
        self.assertTrue(group3 not in user.groups.all())

    def test_login_teams_staff_not_defined(self):
        delattr(settings, 'OPENID_LAUNCHPAD_STAFF_TEAMS')
        user = User.objects.create_user('testuser', 'someone@example.com')
        user.is_staff = True
        user.save()
        self.assertTrue(user.is_staff)

        user = self.get_openid_authed_user_with_teams(user, 'teamname,some-other-team')
        self.assertTrue(user.is_staff)

    def test_login_teams_staff_assignment(self):
        settings.OPENID_LAUNCHPAD_STAFF_TEAMS = ('teamname',)
        user = User.objects.create_user('testuser', 'someone@example.com')
        user.is_staff = False
        user.save()
        self.assertFalse(user.is_staff)

        user = self.get_openid_authed_user_with_teams(user, 'teamname,some-other-team')
        self.assertTrue(user.is_staff)

    def test_login_teams_staff_unassignment(self):
        settings.OPENID_LAUNCHPAD_STAFF_TEAMS = ('different-teamname',)
        user = User.objects.create_user('testuser', 'someone@example.com')
        user.is_staff = True
        user.save()
        self.assertTrue(user.is_staff)

        user = self.get_openid_authed_user_with_teams(user, 'teamname,some-other-team')
        self.assertFalse(user.is_staff)

    def get_openid_authed_user_with_teams(self, user, teams_str):
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()

        # Posting in an identity URL begins the authentication request:
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity'})

        # Complete the request
        openid_request = self.provider.parseFormPost(response.content)
        openid_response = openid_request.answer(True)
        teams_request = teams.TeamsRequest.fromOpenIDRequest(openid_request)
        teams_response = teams.TeamsResponse.extractResponse(
            teams_request, teams_str)
        openid_response.addExtension(teams_response)
        response = self.complete(openid_response)
        return User.objects.get(username=user.username)

    def test_login_complete_signals_login(self):
        # An oauth_login_complete signal is emitted including the
        # request and sreg_response.
        user = User.objects.create_user('someuser', 'someone@example.com')
        useropenid = UserOpenID(
            user=user,
            claimed_id='http://example.com/identity',
            display_id='http://example.com/identity')
        useropenid.save()
        response = self.client.post('/openid/login/',
            {'openid_identifier': 'http://example.com/identity'})
        openid_request = self.provider.parseFormPost(response.content)
        openid_response = openid_request.answer(True)
        # Use a closure to test whether the signal handler was called.
        self.signal_handler_called = False
        def login_callback(sender, **kwargs):
            self.assertTrue(isinstance(
                kwargs.get('request', None), HttpRequest))
            self.assertTrue(isinstance(
                kwargs.get('openid_response', None), SuccessResponse))
            self.signal_handler_called = True
        openid_login_complete.connect(login_callback)

        response = self.complete(openid_response)

        self.assertTrue(self.signal_handler_called)
        openid_login_complete.disconnect(login_callback)


class HelperFunctionsTest(TestCase):
    def test_sanitise_redirect_url(self):
        settings.ALLOWED_EXTERNAL_OPENID_REDIRECT_DOMAINS = [
            "example.com", "example.org"]
        # list of URLs and whether they should be passed or not
        urls = [
            ("http://example.com", True),
            ("http://example.org/", True),
            ("http://example.org/foo/bar", True),
            ("http://example.org/foo/bar?baz=quux", True),
            ("http://example.org:9999/foo/bar?baz=quux", True),
            ("http://www.example.org/", False),
            ("http://example.net/foo/bar?baz=quux", False),
            ("/somewhere/local", True),
            ("/somewhere/local?url=http://fail.com/bar", True),
            # An empty path, as seen when no "next" parameter is passed.
            ("", False),
            ("/path with spaces", False),
        ]
        for url, returns_self in urls:
            sanitised = sanitise_redirect_url(url)
            if returns_self:
                self.assertEqual(url, sanitised)
            else:
                self.assertEqual(settings.LOGIN_REDIRECT_URL, sanitised)

def suite():
    return unittest.TestLoader().loadTestsFromName(__name__)
