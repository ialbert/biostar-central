
from django.shortcuts import render, redirect, reverse
from django.contrib import messages
from . import forms, auth
from .models import Post, Vote, Message
from .decorators import object_exists
from django.contrib.auth.decorators import login_required



def list_view(request, template="forum/post_list.html", extra_context={}, topic=None,
              extra_proc=lambda x:x):
    "List view for posts and messages"

    topic = topic or request.GET.get("topic", 'latest')

    if request.user.is_anonymous and topic != "latest":
        messages.error(request, f"You must be logged in to perform action.")
        topic = "latest"

    objs = auth.list_by_topic(request=request, topic=topic).order_by("-pk")

    # Apply extra protocols to queryset (updates, etc)
    extra_proc(objs)

    context = dict(objs=objs)
    context.update(extra_context)

    return render(request, template_name=template, context=context)


@login_required
def message_list(request):

    active = request.GET.get("q", "inbox")

    amap = dict(inbox="active", outbox="active",unread="active")
    active = active if (active in amap) else "inbox"

    context = {active:amap[active], "not_outbox":active != "outbox"}

    # Change messages in list to "seen",
    # still unread until message_view is visited.
    update_seen = lambda query_set: Message.objects.filter(pk__in=query_set).update(seen=True)

    return list_view(request, template="forum/message_list.html",
                     topic=active, extra_context=context, extra_proc=update_seen)



@object_exists(klass=Message)
def message_view(request, uid):

    return



@object_exists(klass=Post)
@login_required
def update_vote(request, uid):

    # Post to upvote/bookmark
    post = Post.objects.filter(uid=uid).first()
    user = request.user
    vmap = {"upvote": Vote.UP, "bookmark": Vote.BOOKMARK}

    vote_type = vmap.get(request.GET.get("type"), Vote.EMPTY)

    vote = Vote.objects.filter(post=post, author=user, type=vote_type).first()

    if vote:
        # Change vote to empty if clicked twice
        auth.create_vote(update=True, author=user, post=post, vote_type=vote.type,
                         updated_type=Vote.EMPTY)
    elif not vote:
        auth.create_vote(author=user, post=post, vote_type=vote_type)

    return redirect(reverse("post_view", kwargs=dict(uid=post.uid)))



@object_exists(klass=Post)
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


@object_exists(klass=Post)
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




@object_exists(klass=Post)
@login_required
def edit_post(request, uid):
    "Edit an existing post"

    return







