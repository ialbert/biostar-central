

from . import forms, auth

from django.shortcuts import render, redirect, reverse




def post_list(request):
    "List view for posts"


    # show


    return





def post_view(request, uid):
    "Return view for specific post"

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

    return render(request, "post_edit.html", context=context)




def new_answer(request, puid):
    "Provide a short form to make a new answer to a post"


    return





def edit_post(request, uid):
    "Edit an existing post"

    return







