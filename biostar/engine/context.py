from biostar import VERSION
from django.conf import settings


def engine(request):
    '''
    Additional context applied to each request.
    Note: This function is critically important!
    The site will not load up without it.
    '''

    forum_enabaled = settings.ONLY_FORUM_URLS or settings.ENABLE_FORUM
    allow_self_moderate = settings.ALLOW_SELF_MODERATE

    params = dict(user=request.user, VERSION=VERSION, request=request, forum_enabaled=forum_enabaled,
                  allow_self_moderate=allow_self_moderate)

    return params
