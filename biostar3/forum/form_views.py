__author__ = 'ialbert'

from django import forms
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse, reverse_lazy
from django.core.exceptions import ValidationError
from django.contrib.auth.decorators import login_required
from .models import Post
from . import auth
from django.shortcuts import render, redirect

import logging
from django.contrib.auth import get_user_model

logger = logging.getLogger('biostar')

# Get custom user model.
User = get_user_model()

def title_validator(text):
    "Validates form input for tags"
    text = text.strip()
    MIN_LEN, MIN_WORDS = 15, 3
    if not text:
        raise ValidationError('Please enter a title')

    if len(text) < MIN_LEN:
        raise ValidationError('Title is too short! Needs to have at least %s characters' % MIN_LEN)

    words = text.split()
    if len(words) < MIN_WORDS:
        raise ValidationError('Title too simple! Needs more than %s words please.' % MIN_WORDS)


def tag_validator(text):
    MAX_TAGS = 10
    parts = auth.tag_split(text)
    if len(parts) > MAX_TAGS:
        raise ValidationError('Too many tags! Have no more than %s tags please.' % MAX_TAGS)

    if len(parts) < 1:
        raise ValidationError('Please enter at least one tag!')


class ContentForm(forms.Form):
    """
    Edit or create content: answers, comments
    """
    # The is_toplevel field is used to distinguish between subclasses inside templates
    is_toplevel = False
    min_lenght = 5
    content = forms.CharField(widget=forms.Textarea,
                              min_length=min_lenght, max_length=settings.MAX_POST_SIZE,
                              initial="", required=True)


class PostForm(ContentForm):
    """
    Edit or create top level posts: question, news, forum posts,
    """
    is_toplevel = True
    min_lenght = 50

    title = forms.CharField(widget=forms.TextInput, initial='', max_length=200,
                            validators=[title_validator])
    tags = forms.CharField(max_length=100, initial='', validators=[tag_validator])
    type = forms.TypedChoiceField(coerce=int, choices=[
        (Post.QUESTION, "Question"), (Post.NEWS, "News"), (Post.FORUM, "Forum"), (Post.JOB, "Job Ad"),
    ])


def get_post(request, user, pk, edit_access_required=True):
    """
    Authenticates access to a post.
    """
    post = Post.objects.filter(pk=pk).select_related("group", "group__groupinfo").first()

    if not post:
        messages.error(request, "Post does not exist. Perhaps it has been deleted.")
        raise auth.AccessDenied()

    if not auth.read_access_post(user=user, post=post):
        messages.error(request, "This post may not be accessed by this user!")
        raise auth.AccessDenied()

    if edit_access_required and not auth.write_access_post(user, post):
        messages.error(request, "This post may not be edited by this user!")
        raise auth.AccessDenied()

    return post


@login_required
def create_node(request, parent_id=None, post_type=None):
    """
    This view creates nodes. Is not called directly from the web only through
    other functions that prefill parameters.
    """
    user = request.user
    group = request.group
    template_name = "edit_post.html"
    redirect_home = redirect(reverse_lazy("home"))

    # No post type means a top level post to be created.
    if post_type is None:
        form_class = PostForm
        action = reverse("new_post")
    else:
        form_class = ContentForm
        if post_type == Post.ANSWER:
            action = reverse("new_answer", kwargs=dict(parent_id=parent_id))
        else:
            action = reverse("new_comment",  kwargs=dict(parent_id=parent_id))

    if request.method == "GET":
        # This will render the initial form for the user.
        if parent_id is not None:
            try:
                # Need to make sure that the parent post is readable to the user.
                parent = get_post(request=request, user=user, pk=parent_id, edit_access_required=False)
            except auth.AccessDenied:
                return redirect_home

        form = form_class()
        context = dict(form=form, action=action)
        return render(request, template_name, context)

    if request.method == "POST":
        # Data is being submitted
        form = form_class(request.POST)

        if parent_id is not None:
            # Need to make sure that the parent post is readable to the user.
            try:
                parent = get_post(request=request, user=user, pk=parent_id, edit_access_required=False)
            except auth.AccessDenied:
                return redirect_home

        if not form.is_valid():
            # Form data came but not valid.
            context = dict(form=form, action=action)
            return render(request, template_name, context)

        # The form is valid create the post based on the form.
        if post_type is None:
            post = auth.create_toplevel_post(user=user, group=group, data=form.cleaned_data)
        else:
            post = auth.create_content_post(data=form.cleaned_data, post_type=post_type, user=user, parent=parent)

        return redirect(post.get_absolute_url())


def create_toplevel_post(request):
    "A new toplevel post"
    return create_node(request=request, parent_id=None, post_type=None)


def create_answer(request, parent_id):
    return create_node(request=request, parent_id=parent_id, post_type=Post.ANSWER)


def create_comment(request, parent_id):
    return create_node(request=request, parent_id=parent_id, post_type=Post.COMMENT)



@login_required
def edit_post(request, post_id):
    """
    This view updates posts.
    """
    user = request.user
    template_name = "edit_post.html"
    action = reverse("edit_post", kwargs=dict(post_id=post_id))
    redirect_home = redirect(reverse_lazy("home"))

    try:
        # Need to make sure that the parent post is readable to the user.
        post = get_post(request=request, user=user, pk=post_id, edit_access_required=True)
    except auth.AccessDenied:
        return redirect_home

    # Different forms chosen based on post type.
    if post.is_toplevel:
        form_class = PostForm
        tags = ", ".join(post.tags.names())
        initial = dict(content=post.content, title=post.title, tags=tags, type=post.type)
    else:
        form_class = ContentForm
        initial = dict(content=post.content)


    if request.method == "GET":
        form = form_class(initial=initial)
        context = dict(form=form, action=action)
        return render(request, template_name, context)

    if request.method == "POST":
        form = form_class(request.POST)

        # Invalid form, bail out with error messaged.
        if not form.is_valid():
            context = dict(form=form, action=action)
            return render(request, template_name, context)

        # The data is valid update the post and return the view to it.
        get = lambda word: form.cleaned_data.get(word, '')

        post.content = get('content')
        post.lastedit_user = user
        post.lastedit_date = auth.now()

        # Extra information to be saved for toplevel posts.
        if post.is_toplevel:
            post.title = get('title')
            post.type = get('type')
            tags = get('tags')
            tags = auth.tag_split(tags)
            # Must explicitly set the new tags.
            post.tags.set(*tags)

        post.save()

    return redirect(post.get_absolute_url())
