from django.shortcuts import render
from django.db.models import Max, Count
from django.core.paginator import Paginator
from biostar.planet.models import Blog, BlogPost
from django.conf import settings
from biostar.utils.decorators import reset_count


@reset_count(key="planet_count")
def blog_list(request):

    page = request.GET.get("page", 1)
    blogposts = BlogPost.objects.select_related("blog").order_by("-creation_date")

    blogs = Blog.objects.annotate(updated_date=Max("blogpost__creation_date"))
    blogs = blogs.annotate(count=Count("blogpost__id"))
    blogs = blogs.order_by("-updated_date", "-list_order")[:100]

    # .distinct() on textfield not allowed on sql database backend.
    # .distinct() + .annotate() not implemented in postgres unless ordered by text field.
    # seen = set()
    #
    # def distinct(blg):
    #     # Return True if feed is already seen
    #     if blg.feed in seen:
    #         return False
    #     seen.update([blg.feed])
    #     return True
    #
    # blogs = [b for b in blogs if distinct(b)]

    blogposts = Paginator(blogposts, per_page=settings.BLOGS_PER_PAGE)
    blogposts = blogposts.get_page(page)

    context = dict(blogposts=blogposts, tab='planet', blogs=blogs)
    return render(request, 'planet/blog_list.html', context)

