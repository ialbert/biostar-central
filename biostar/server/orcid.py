import json

from requests_oauthlib import OAuth2Session
from allauth.socialaccount.models import SocialApp


def import_orcid_profile_data(**kwargs):
    provider = kwargs['sociallogin'].account.provider

    if not 'orcid' in provider.lower():
        return

    social_token = kwargs['sociallogin'].token
    uid = kwargs['sociallogin'].account.uid
    data = _run_profile_query(provider, social_token.token, uid)

    print(json.dumps(data, indent=4))


def _run_profile_query(provider, token, uid):
    oauth_sex = _create_orcid_oauth_sex(provider, token)
    resource_url = 'https://api.sandbox.orcid.org/v1.1/{}/orcid-profile'.format(uid)
    r = oauth_sex.get(resource_url, headers={'accept': 'application/orcid+json'})
    return r.json()


def _create_orcid_oauth_sex(provider, token):
    orcid_provider = SocialApp.objects.get(provider=provider)
    sex = OAuth2Session(client_id=orcid_provider.client_id,
                        token={'access_token': token,})
    return sex