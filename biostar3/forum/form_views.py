__author__ = 'ialbert'

from django import forms
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.core.exceptions import ValidationError
from braces.views import LoginRequiredMixin
from django.views.generic import FormView
from django.contrib.auth.decorators import login_required
from django.views.decorators.csrf import csrf_exempt
from django.core.exceptions import ObjectDoesNotExist
from .models import Post
from django.http import HttpResponseRedirect
from . import auth

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

    content = forms.CharField(widget=forms.Textarea,
                              min_length=settings.MIN_POST_SIZE, max_length=settings.MAX_POST_SIZE,
                              initial="", required=True)


class PostForm(ContentForm):
    """
    Edit or create top level posts: question, news, forum posts,
    """
    title = forms.CharField(widget=forms.TextInput, initial='', max_length=200,
                            validators=[title_validator])
    tags = forms.CharField(max_length=100, initial='', validators=[tag_validator])
    type = forms.TypedChoiceField(coerce=int, choices=[
        (Post.QUESTION, "Question"), (Post.NEWS, "News"), (Post.FORUM, "Forum"), (Post.JOB, "Job Ad"),
    ])


def redirect(name, **kwargs):
    return HttpResponseRedirect(reverse(name, kwargs=kwargs))

class BaseNode(LoginRequiredMixin, FormView):
    form_class = ContentForm
    template_name = "edit_node.html"
    edit_access_required = True

    def get_post(self, user):
        """
        Will authenticate user access to
        """
        pk = self.kwargs.get('pk')
        post = Post.objects.filter(pk=pk).select_related("group", "group__groupinfo").first()

        if not post:
            messages.error(self.request, "Post does not exist. Perhaps it has been deleted.")
            raise auth.AccessDenied()

        if not auth.read_access_post(user=user, post=post):
            messages.error(self.request, "This post may not be accessed by this user!")
            raise auth.AccessDenied()

        if self.edit_access_required and not auth.write_access_post(user, post):
            messages.error(self.request, "This post may not be edited by this user!")
            raise auth.AccessDenied()

        return post

    def action(self, post, user, form):
        raise NotImplementedError()

    def get(self, *args, **kwargs):
        try:
            self.get_post(user=self.request.user)
        except auth.AccessDenied, exc:
            logger.error(exc)
            return redirect("home")

        return super(BaseNode, self).get(*args, **kwargs)

    def post(self, *args, **kwargs):
        try:
            self.get_post(user=self.request.user)
        except auth.AccessDenied, exc:
            logger.error(exc)
            return redirect("home")
        return super(BaseNode, self).post(*args, **kwargs)

    def form_valid(self, form):
        # Create the new post.
        user = self.request.user
        post = self.get_post(user=user)
        content = form.cleaned_data['content']
        try:
            # Apply the action using the content and post
            self.post = self.action(post=post, user=user, form=form)
        except Exception, exc:
            logger.error(exc)
            messages.error(self.request, "Unable to create the post!")
            return redirect("home")

        return super(BaseNode, self).form_valid(form)

    def get_success_url(self):
        return self.post.get_absolute_url()


class NewNode(BaseNode):
    edit_access_required = False

    def action(self, post, user, form):
        # The incoming post is the parent in this case.
        content = form.cleaned_data.get('content', '')
        obj = Post.objects.create(parent=post, content=content, author=user)
        return obj


class EditNode(BaseNode):

    def get(self, *args, **kwargs):

        try:
            post = self.get_post(user=self.request.user)
            if post.is_toplevel:
                return redirect("edit_post", pk=post.id)

        except auth.AccessDenied, exc:
            logger.error(exc)
            return redirect("home")

        return super(EditNode, self).get(*args, **kwargs)

    def get_initial(self):
        try:
            post = self.get_post(user=self.request.user)
        except auth.AccessDenied, exc:
            logger.error(exc)
            return redirect("home")
        initial = dict(content=post.content)
        return initial

    def action(self, post, user, form):
        # Update the post
        post.content = form.cleaned_data.get('content', '')
        post.lastedit_user = user
        post.lastedit_date = auth.now()
        post.save()
        return post


class NewPost(BaseNode):
    form_class = PostForm
    template_name = "edit_post.html"
    edit_access_required = False

    def get_post(self, user):
        return None

    def action(self, post, user, form):
        # The incoming post is None in this case.

        title = form.cleaned_data.get('title', '').strip()
        type = form.cleaned_data.get('type', '')
        tags = form.cleaned_data.get('tags', '')
        tags = auth.tag_split(tags)
        content = form.cleaned_data.get('content', '')

        obj = Post.objects.create(content=content, title=title,
                                  author=user, type=type, group=self.request.group)

        # Set the tags on the post
        obj.tags.set(*tags)

        # Self referential ForeignKeys need to be updated explicitly!
        Post.objects.filter(pk=obj.pk).update(root_id=obj.id, parent_id=obj.id)

        return obj

class EditPost(BaseNode):
    form_class = PostForm
    template_name = "edit_post.html"

    def get_initial(self):

        try:
            post = self.get_post(user=self.request.user)
        except auth.AccessDenied, exc:
            logger.error(exc)
            return redirect("home")

        tags = ", ".join(post.tags.names())
        initial = dict(content=post.content, title=post.title, tags=tags, type=post.type)
        return initial

    def action(self, post, user, form):
        get = form.cleaned_data.get
        post.content = get('content', '')
        post.title = get('title', '').strip()
        post.type = get('type', '')
        tags = get('tags', '')
        tags = auth.tag_split(tags)
        post.lastedit_user = user
        post.lastedit_date = auth.now()
        post.save()

        # Set the new tags.
        post.tags.set(*tags)

        return post