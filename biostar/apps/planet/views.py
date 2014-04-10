# Create your views here.
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, View
from .models import Blog, BlogPost

class BlogPostList(ListView):
    template_name = "planet/planet_entries.html"
    paginate_by = 25
    model = BlogPost
    context_object_name = 'blogposts'

    def get_queryset(self):
        query = super(BlogPostList, self).get_queryset()
        return query.select_related("blog").order_by("-creation_date")



