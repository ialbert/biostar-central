

from . import forms


def new_post(request):
    "Make a new post"

    form = forms.PostLongForm()

    if request.method == "POST":
        form = forms.PostLongForm(data=request.POST)
        if form.is_valid():
            # Create a new post by user
            form.save(author=request.user)
            return

    return





def new_answer(request, puid):
    "Make a new answer to a post"


    return





def edit_post(request, uid):
    "Edit a post"

    return







