from django.shortcuts import render_to_response
from django.views import generic
from biostar.apps.accounts.models import User

# Create your views here.
class IndexView(generic.TemplateView):
    name="index"
    template_name = "index.html"

    def get_context_data(self, **kwargs):
        context = super(IndexView, self).get_context_data(**kwargs)
        users = User.objects.all()

        users[0].name = "A"
        users[0].save()

        return context