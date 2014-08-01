import json
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import redirect


def hook_social_account_added(**kwargs):
    """
    Method executed when a new user is created using a social provider or a new social
    provider is connected to an existing user.

    `kwargs` is a dictionary in the following form:
    {
      'signal': <django.dispatch.dispatcher.Signal object at ...>,
      'sociallogin': <allauth.socialaccount.models.SocialLogin object at ...>,
      'request': <WSGIRequest ...>,
      'sender': <class 'allauth.socialaccount.models.SocialLogin'>
    }
    """
    if 'orcid' in kwargs['sociallogin'].account.provider.lower():
        ask_to_import_orcid_profile(kwargs['request'])


def ask_to_import_orcid_profile(request):
    """
    Ask the user if he want to import his ORCID profile data.

    Parameters:
    request - a `WSGIRequest`;
    """
    messages.info(request, "Do you want <a href='{}'>one</a>?".format(reverse('orcid-import')))


@login_required
def import_bio(request):
    """
    Import bio from ORCID for the current user and store it in ........
    """
    # The ORCID bio for the current user is stored in:
    # `allauth.socialaccount.models.SocialAccount.extra_data`
    # related to the current user.
    user = request.user
    social_account = user.socialaccount_set.get(provider__icontains='orcid')
    data = social_account.extra_data
    works = extract_from_dict(data, ['orcid-profile', 'orcid-activities', 'orcid-works',
                                     'orcid-work'])
    # Add works (list of papers).
    txt = '<b>Works</b><br /><ul>'
    for work in works[:3]:
        title = extract_from_dict(work, ['work-title', 'title', 'value'])
        year = extract_from_dict(work, ['publication-date', 'year', 'value'])
        month = extract_from_dict(work, ['publication-date', 'month', 'value'])
        day = extract_from_dict(work, ['publication-date', 'day', 'value'])
        txt += '<li><i>{}</i>{}{}{}</li>'.format(
            title,
            ', ' + year if year and month else '',
            '/' + month if year and month else '',
            '/' + day if year and month and day else ''
        )
    txt += '</ul>'

    # Add bio.
    bio = extract_from_dict(data, ['orcid-profile', 'orcid-bio', 'biography', 'value'])
    if bio:
        bio += '<br />'
    txt = bio + txt

    user.profile.info = txt
    user.profile.save()
    return redirect(reverse('user-details', kwargs={'pk': user.pk}))


def extract_from_dict(data, path):
    """
    Navigate `data`, a multidimensional array (list or dictionary), and return the object
    at `path`.
    """
    value = data
    try:
        for key in path:
            value = value[key]
        return value
    except Exception:
        return ''