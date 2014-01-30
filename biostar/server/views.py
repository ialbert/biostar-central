from django.shortcuts import render_to_response
from django.views import generic
from biostar.apps.people.models import User


class PageBase(generic.TemplateView):
    """
    All pages will inherit the attributes of this class.
    """
    page_title = "Bioinformatics Answers on Biostars"

    def get_context_data(self, **kwargs):
        context = super(PageBase, self).get_context_data(**kwargs)
        context['page_title'] = self.page_title
        context['user'] = self.request.user
        return context


class IndexView(PageBase):
    page_title = "Welcome to Biostars"
    template_name = "index.html"


class UserView(PageBase):
    page_title = "Users on Biostars"
    template_name = "users.html"