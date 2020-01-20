from django.shortcuts import render
from django.db.models import Max, Count
from django.core.paginator import Paginator
from biostar.planet.models import Blog, BlogPost
from django.conf import settings


def blog_list(request):

    page = request.GET.get("page", 1)
    blogposts = BlogPost.objects.select_related("blog").order_by("-creation_date")
    blogs = Blog.objects.all().annotate(updated_date=Max("blogpost__creation_date"),
                                        count=Count("blogpost__id")).order_by("-updated_date", "-list_order")
    blogposts = Paginator(blogposts, per_page=settings.BLOGS_PER_PAGE)
    blogposts = blogposts.get_page(page)

    context = dict(blogposts=blogposts, tab='planet', blogs=blogs)
    return render(request, 'planet/blog_list.html', context)

