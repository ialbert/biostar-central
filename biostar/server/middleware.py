__author__ = 'ialbert'

class Visit(object):
    """
    Sets visit specific parameters on objects.
    """

    def process_request(self, request):

        user = request.user
        if not user.is_authenticated():
            # This attribute is required inside templates.
            user.is_moderator = False

