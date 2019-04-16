from biostar.forum import tasks


def forum_middleware(get_response):

    def middleware(request):

        user = request.user

        # Handle tasks async
        if tasks.HAS_UWSGI and user.is_authenticated:
            tasks.create_user_awards(user_id=user.id)

        response = get_response(request)
        # Can process response here after its been handled by the view

        return response

    return middleware
