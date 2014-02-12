__author__ = 'ialbert'
from braces.views import JSONResponseMixin
from biostar.apps.posts.models import Post
from biostar.apps.users.models import User
from django.views.generic import View


class VoteSubmit(JSONResponseMixin, View):
    json_dumps_kwargs = {'indent': 2}

    def get(self, request, *args, **kwargs):

        try:
            context_dict = {
                'title': "Hello world!",
            }
        except Exception, exc:
            context_dict = {
                'error': "Error"
            }

        return self.render_json_response(context_dict)
