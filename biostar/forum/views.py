
from django.shortcuts import render, redirect, reverse
from django.contrib import messages
from . import forms, auth
from .models import Post, Vote
from .decorators import object_exists

LATEST = "Latest"


def post_list(request):
    "List view for posts"

    posts = auth.posts_by_topic(request=request, topic=LATEST)

    context = dict(posts=posts)

    return render(request, "forum/post_list.html", context=context)

@object_exists
def post_view(request, uid):
    "Return a detailed view for specific post"

    user = request.user

    obj = Post.objects.filter(uid=uid).first()

    # Raise 404 if a deleted post is viewed by an anonymous user
    if (obj.status == Post.DELETED) and not user.profile.is_moderator:
        return redirect(reverse("post_list"))

    # Update the post views.
    Post.update_post_views(obj, request=request)

    # Adds the permissions
    obj = auth.post_permissions(request=request, post=obj)

    # This will be piggybacked on the main object.
    #post.sub = Subscription.get_sub(post=obj, user=user)

    # Return parent view if not at top level.
    if not obj.is_toplevel:
        return redirect(reverse("post_view", kwargs=dict(uid=obj.uid)))

    # Populate the object to build a tree that contains all posts in the thread.
    # Answers sorted before comments.
    thread = [auth.post_permissions(request=request, post=p)
              for p in Post.objects.get_thread(obj, user)]

    # Build tree and gather votes.
    tree = auth.build_tree(thread=thread)
    votes = auth.get_votes(user=user, thread=thread)

    # Shortcuts to each storage.
    bookmarks = votes[Vote.BOOKMARK]
    upvotes = votes[Vote.UP]

    def decorate(post):
        # Can the current user accept answers
        post.has_bookmark = post.id in bookmarks
        post.has_upvote = post.id in upvotes
        post.can_accept = obj.author == user or post.has_accepted

    # Add attributes by mutating the objects
    map(decorate, thread + [obj])

    # Additional attributes used during rendering
    obj.tree = tree
    obj.answers = [p for p in thread if p.type == Post.ANSWER]

    context = dict(post=obj)
    return render(request, "forum/post_view.html", context=context)




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



@object_exists
def new_answer(request, puid):
    "Provide a short form to make a new answer to a post"


    return




@object_exists
def edit_post(request, uid):
    "Edit an existing post"

    return







