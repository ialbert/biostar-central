from django.shortcuts import render_to_response
from django.views import generic

# Create your views here.
class IndexView(generic.TemplateView):
    name="index"
    template_name = "index.html"