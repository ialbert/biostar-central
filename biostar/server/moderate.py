"""
Moderator views
"""
from django.conf import settings
from django.views.generic import DetailView, ListView, TemplateView, RedirectView, View
from biostar.apps.posts.auth import post_permissions

class ModeratorPanel(TemplateView):
    template_name = "moderator-panel.html"

    def post(self, *args, **kwargs):
        return