from django.shortcuts import render_to_response
from django.views import generic
from biostar.apps.people.models import User

# Create your views here.
class IndexView(generic.TemplateView):
    name = "index"
    page_title = "Bioinformatics Answers on Biostars"
    template_name = "index.html"

    def get_context_data(self, **kwargs):
        context = super(IndexView, self).get_context_data(**kwargs)
        context['page_title'] = self.page_title
        return context