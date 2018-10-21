from biostar import VERSION
from django.conf import settings
from biostar.engine.const import *


def engine(request):
    '''
    Additional context applied to each request.
    Note: This function is critically important!
    The site will not load up without it.
    '''
    enable_forum = settings.ENABLE_FORUM or settings.ONLY_FORUM_URLS

    params = dict(user=request.user, VERSION=VERSION, request=request,
                data_board=DATA_CLIPBOARD, recipe_board=RECIPE_CLIPBOARD,
                results_board=RESULTS_CLIPBOARD, enable_forum=enable_forum,
                only_enable_forum=settings.ONLY_FORUM_URLS
    )

    return params
