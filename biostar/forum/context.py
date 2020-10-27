from biostar import VERSION
from django.conf import settings
def forum(request):
    '''
    Additional context applied to each request.
    Note: This function is critically important!
    The site will not load up without it.
    '''

    #print(request.session.get('res'), 'res', request.COOKIES)
    res = request.COOKIES.get('resolution', 'x')
    width, height = res.split('x')

    params = dict(user=request.user, width=width, height=height,
                  VERSION=VERSION, request=request, site_name=settings.SITE_NAME,
                  site_domain=settings.SITE_DOMAIN, google_tracker=settings.GOOGLE_TRACKER)
    return params
