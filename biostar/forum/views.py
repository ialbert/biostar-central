
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
@login_required
def update_vote(request, uid):

    # Post to upvote/bookmark
    post = Post.objects.filter(uid=uid).first()
    user = request.user

    vote_type = request.GET.get("type", Vote.EMPTY)

    vote = Vote.objects.filter(post=post, author=user).first()

    # Pressing vote button multiple times toggles
    cond = (vote.type == vote_type and vote_type)
    vote_type = Vote.EMPTY if cond else vote_type

    if vote.exists():
        # Update an existing vote
        auth.create_vote(update=True,author=user, post=post, vote_type=vote_type)
    else:
        auth.create_vote(author=user, post=post, vote_type=vote_type)

    return redirect(reverse("post_view", kwargs=dict(uid=post.uid)))



@object_exists
def post_view(request, uid):
    "Return a detailed view for specific post"

    # Form used for answers
    form = forms.PostShortForm()

    # Get the parents info
    obj = Post.objects.filter(uid=uid).first()

    # Return root view if not at top level.
    obj = obj if obj.is_toplevel else obj.root

    # Update the post views.
    Post.update_post_views(obj, request=request)

    if request.method == "POST":
        form = forms.PostShortForm(data=request.POST)
        if form.is_valid():
            form.save(parent=obj.parent, author=request.user)
            return redirect(reverse("post_view", kwargs=dict(uid=obj.root.uid)))

    # Adds the permissions
    obj = auth.post_permissions(request=request, post=obj)

    # Populate the object to build a tree that contains all posts in the thread.
    # Answers are added here as well.
    obj = auth.build_obj_tree(request=request, obj=obj)

    context = dict(post=obj, form=form)
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


@object_exists
@login_required
def subs_action(request, uid):

    # Post actions are being taken on
    post = Post.objects.filter(uid=uid).first()
    user = request.user

    if request.method == "POST" and user.is_authenticated:
        form = forms.SubsForm(data=request.POST, post=post, user=user)

        if form.is_valid():
            sub = form.save()
            msg = f"Updated Subscription to : {sub.get_type_display()}"
            messages.success(request, msg)

    return redirect(reverse("post_view", kwargs=dict(uid=post.uid)))


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







