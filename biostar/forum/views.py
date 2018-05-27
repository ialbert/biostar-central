

from . import forms, auth

from django.shortcuts import render, redirect, reverse

LATEST = "Latest"


def post_list(request):
    "List view for posts"

    posts = auth.posts_by_topic(request=request, topic=LATEST)

    context = dict(posts=posts)

    return render(request, "forum/post_list.html", context=context)


def post_view(request, uid):
    "Return a detailed view for specific post"

    return




def new_post(request):
    "Make a new post"

    form = forms.PostLongForm()

    if request.method == "POST":
        form = forms.PostLongForm(data=request.POST)
        if form.is_valid():
            # Create a new post by user
            post = form.save(author=request.user)
            return redirect(reverse("post_view", kwargs=dict(uid=post.uid)))

    context = dict(form=form)

    return render(request, "forum/post_edit.html", context=context)




def new_answer(request, puid):
    "Provide a short form to make a new answer to a post"


    return





def edit_post(request, uid):
    "Edit an existing post"

    return







