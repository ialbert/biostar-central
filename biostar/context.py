from biostar import VERSION
from django.conf import settings

def main(request):
    '''
    Additional context applied to each request.
    Note: This function is critically important!
    The site will not load up without it.
    '''

    params = dict(user=request.user, VERSION=VERSION, site_name=settings.SITE_NAME, request=request)

    return params
