from biostar import VERSION
from django.conf import settings

def engine(request):
    '''
    Additional context applied to each request.
    Note: This function is critically important!
    The site will not load up without it.
    '''
    params = dict(user=request.user, VERSION=VERSION, SITE_HEADER=settings.SITE_HEADER)

    return params
