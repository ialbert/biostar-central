from django.shortcuts import render
from django.db.models import Max, Count
from django.core.paginator import Paginator
from biostar.planet.models import Blog, BlogPost
from django.conf import settings
from biostar.utils.decorators import reset_count
from django.shortcuts import reverse
from django.shortcuts import render, redirect
from django.contrib import messages
from django.http import Http404
from biostar.forum import util

@reset_count(key="planet_count")
def blog_list(request):

    page = request.GET.get("page", 1)
    blogposts = BlogPost.objects.select_related("blog").order_by("-rank", "-creation_date")

    blogs = Blog.objects.annotate(updated_date=Max("blogpost__creation_date"))
    blogs = blogs.annotate(count=Count("blogpost__id"))
    blogs = blogs.order_by("-updated_date", "-list_order")[:100]

    blogposts = Paginator(blogposts, per_page=settings.BLOGS_PER_PAGE)
    blogposts = blogposts.get_page(page)

    context = dict(blogposts=blogposts, tab='planet', blogs=blogs)
    return render(request, 'planet/blog_list.html', context)

def blog_bump(request, id):

    post = BlogPost.objects.filter(id=id).update(rank=util.now())

    return redirect(to=reverse('blog_list'))


def blog_view(request, id):

    post = BlogPost.objects.filter(id=id).first()

    if not post:
        msg = f"Invalid blog post id: {id}"
        raise Http404(msg)

    # Redirect to remote
    if post.blog.remote:
        return redirect(post.link)

    context = dict(post=post)

    return render(request, 'planet/blog_view.html', context)

