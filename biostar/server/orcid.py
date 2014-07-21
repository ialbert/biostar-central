"""
ORCID

Info: http://support.orcid.org/knowledgebase/articles/180285-introduction-to-the-orcid-api
Register a client app: http://orcid.org/content/register-client-application

Production API
==============
Urls:
    - Public: http://pub.orcid.org
    - Member: https://api.orcid.org

Sandbox (test) API
==================
Info: http://support.orcid.org/knowledgebase/articles/166623-about-the-orcid-sandbox
Urls:
    - Public: http://pub.api.sandbox.orcid.org
    - Member: http://api.sandbox.orcid.org

Create some members here: https://sandbox.orcid.org/register
test@test123.com
http://sandbox.orcid.org/0000-0001-6796-198X

API endpoints and resources
===========================
http://support.orcid.org/knowledgebase/articles/116874-orcid-api-guide

Access Tokens and Scopes
========================
http://support.orcid.org/knowledgebase/articles/117203-obtaining-access-tokens
http://support.orcid.org/knowledgebase/articles/120162-orcid-scopes

OAuth
=====
- Help: http://support.orcid.org/knowledgebase/articles/154641-introduction-to-the-api-with-google-s-oauth-playgr
- Step1: http://support.orcid.org/knowledgebase/articles/120107-get-oauth-authorize
- Step2: http://support.orcid.org/knowledgebase/articles/119985
- Example: http://support.orcid.org/knowledgebase/articles/179969-methods-to-generate-an-access-token-for-testing
- Playground: https://developers.google.com/oauthplayground/#step3&url=https%3A//api.sandbox.orcid.org/&content_type=application/json&http_method=GET&useDefaultOauthCred=unchecked&oauthEndpointSelect=Custom&oauthAuthEndpointValue=https%3A//sandbox.orcid.org/oauth/authorize&oauthTokenEndpointValue=https%3A//api.sandbox.orcid.org/oauth/token&includeCredentials=unchecked&accessTokenType=bearer&autoRefreshToken=unchecked&accessType=offline&forceAprovalPrompt=checked&response_type=code

Examples
========
Public API (no access token required): http://pub.sandbox.orcid.org/v1.1/0000-0001-6796-198X/orcid-bio
"""
from django.http import HttpResponse
from django.shortcuts import redirect
from django.conf import settings

from requests_oauthlib import OAuth2Session
import requests
import json


def test(request):
    html = "<html><a href='authorize/'>go</a></html>"
    return HttpResponse(html)


def test_public_api():
    r = requests.get('http://pub.sandbox.orcid.org/v1.1/0000-0001-6796-198X/orcid-bio',
                     headers={'accept': 'application/orcid+json'})

    print(json.dumps(r.json(), indent=4))

    print("{}Given names: {}".format(
        "\n\n**\n",
        r.json()['orcid-profile']['orcid-bio']['personal-details']['given-names']['value']))

    print("ORCID identifier: {}".format(
        r.json()['orcid-profile']['orcid-identifier']['path']))


def step1(request):
    oauth_sex = _create_oauth_sex()
    authorization_url, csrf_code = _get_authorization_url(oauth_sex)

    # Store csrf_code to prevent CSRF (OAuth1 has no such code)
    if csrf_code:
        request.session['csrf_code'] = csrf_code

    print("Auth url: {}".format(authorization_url))
    print("Code: {}".format(csrf_code))
    return redirect(authorization_url)


def _create_oauth_sex(token=None):
    if not token:
        sex = OAuth2Session(
            client_id=settings.ORCID_CLIENT_ID,
            scope=u'/orcid-profile/read-limited',  # Note: must be unicode.
            redirect_uri='http://127.0.0.1:8000/orcid-test/callback',
        )
    else:
        sex = OAuth2Session(
            client_id=settings.ORCID_CLIENT_ID,
            token=token)
    return sex


def _get_authorization_url(oauth_sex):
    authorization_url, csrf_code = oauth_sex.authorization_url(
        'https://sandbox.orcid.org/oauth/authorize',
    )
    return authorization_url, csrf_code


def step2(request):
    _prevent_csrf(request)

    auth_code = request.GET.get('code', None)
    sex = _create_oauth_sex()
    token_set = _fetch_token(sex, auth_code)
    # `token_set` has the following form:
    # {
    #   u'access_token': u'...',
    #   u'expires_in': 631138518,
    #   u'expires_at': 2037116888.678718,
    #   u'token_type': u'bearer',
    #   u'orcid': u'0000-0001-6796-198X',
    #   u'scope': [u'/orcid-profile/read-limited'],
    #   u'refresh_token': u'...'
    # }

    print("Token set: {}".format(token_set))

    txt = _test_query(token_set)
    return HttpResponse("<html><pre>{}</pre></html>".format(json.dumps(txt, indent=4)))


def _prevent_csrf(request):
    """Check csrf_code to prevent CSRF"""

    # The `state` GET argument received from the provider in step2 must match the `csrf_code`
    # stored in the user session during step1

    csrf_code = request.GET.get('state', 'default')

    if csrf_code != request.session.pop('csrf_code', None):
        # Raise an exception
        # This will generate a HTTP 500 Server error page
        msg = ("The state parameter received in step2 doesn't match the one in step1. "
               "This can be a security issue.")
        raise Exception(msg)


def _fetch_token(oauth_sex, auth_code):
    return oauth_sex.fetch_token(
        'https://api.sandbox.orcid.org/oauth/token',
        client_secret=settings.ORCID_CLIENT_SECRET,
        code=auth_code
    )


def _test_query(token_set):
    oauth_sex = _create_oauth_sex(token=token_set)
    resource_url = 'https://api.sandbox.orcid.org/v1.1/{}/orcid-profile'.format(token_set['orcid'])
    r = oauth_sex.get(resource_url,
                      headers={'accept': 'application/orcid+json'})
    return r.json()