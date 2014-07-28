import json


def import_orcid_profile_data(**kwargs):
    # The ORCID bio for a user was stored in the `extra_data` of his `SocialAccount` when he added
    # ORCID as social login.
    data = kwargs['sociallogin'].account.extra_data
    print(json.dumps(data, indent=4))