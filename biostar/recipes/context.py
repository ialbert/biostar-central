from biostar import VERSION


def engine(request):
    '''
    Additional context applied to each request.
    Note: This function is critically important!
    The site will not load up without it.
    '''

    params = dict(user=request.user, VERSION=VERSION, request=request)

    return params
