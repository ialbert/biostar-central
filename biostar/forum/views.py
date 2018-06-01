
from django.shortcuts import render, redirect, reverse
from django.contrib import messages
from . import forms, auth
from .models import Post, Vote
from .decorators import object_exists
from django.contrib.auth.decorators import login_required




def post_list(request):
    "List view for posts"

    topic = request.GET.get("topic", 'latest')

    if request.user.is_anonymous and topic != "latest":
        messages.error(request, f"You must be logged in to perform action.")
        topic = "latest"

    posts = auth.posts_by_topic(request=request, topic=topic).order_by("-pk")

    context = dict(posts=posts)

    return render(request, "forum/post_list.html", context=context)



@object_exists
def post_view(request, uid):
    "Return a detailed view for specific post"

    user = request.user

    # Form used for answers
    form = forms.PostShortForm()

    # Get the parents info
    obj = Post.objects.filter(uid=uid).first()

    add_comment = request.GET.get("add", '') == 'comment'

    # Raise 404 if a deleted post is viewed by an anonymous user
    if (obj.status == Post.DELETED) and not user.profile.is_moderator:
        return redirect(reverse("post_list"))

    # Return root view if not at top level.
    if not obj.is_toplevel:
        return redirect(reverse("post_view", kwargs=dict(uid=obj.root.uid)))

    # Update the post views.
    Post.update_post_views(obj, request=request)

    if request.method == "POST":
        form = forms.PostShortForm(data=request.POST)
        if form.is_valid():
            form.save(parent=obj.parent, author=request.user)
            return redirect(reverse("post_view", kwargs=dict(uid=obj.root.uid)))

    # Adds the permissions
    obj = auth.post_permissions(request=request, post=obj)

    # This will be piggybacked on the main object.
    #post.sub = Subscription.get_sub(post=obj, user=user)

    # Populate the object to build a tree that contains all posts in the thread.
    # Answers are added here as well.
    obj = auth.build_obj_tree(request=request, obj=obj)

    context = dict(post=obj, form=form, add_comment=add_comment)
    return render(request, "forum/post_view.html", context=context)


@login_required
def post_comment(request, uid):

    # Get the parent post to add comment to
    obj = Post.objects.filter(uid=uid).first()

    # Form used for answers
    form = forms.PostShortForm()

    if request.method == "POST":

        form = forms.PostShortForm(data=request.POST)

        if form.is_valid():
            form = forms.PostShortForm(data=request.POST)
            if form.is_valid():
                form.save(parent=obj, author=request.user, post_type=Post.COMMENT)
            return redirect(reverse("post_view", kwargs=dict(uid=obj.root.uid)))

    context = dict(form=form, post=obj)
    return render(request, "forum/post_comment.html", context=context)



@login_required
def post_create(request):
    "Make a new post"

    form = forms.PostLongForm()

    if request.method == "POST":
        form = forms.PostLongForm(data=request.POST)
        if form.is_valid():
            # Create a new post by user
            post = form.save(author=request.user)
            return redirect(reverse("post_view", kwargs=dict(uid=post.uid)))

    context = dict(form=form)

    return render(request, "forum/post_create.html", context=context)




@object_exists
@login_required
def edit_post(request, uid):
    "Edit an existing post"

    return







