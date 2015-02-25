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
    title = forms.CharField(widget=forms.TextInput, initial='', max_length=200)
    tags = forms.CharField(max_length=100, initial='')
    type = forms.TypedChoiceField(coerce=int, choices=[
        (Post.QUESTION, "Question"), (Post.NEWS, "News"), (Post.FORUM, "Forum"),(Post.JOB, "Job Ad"),
    ])

def redirect(name):
    return HttpResponseRedirect(reverse(name))


class BaseNode(LoginRequiredMixin, FormView):
    form_class = ContentForm
    template_name = "edit_node.html"

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
    def action(self, post, user, form):
        content = form.cleaned_data.get('content', '')
        obj = Post.objects.create(parent=post, content=content, author=user)
        return obj


class EditNode(BaseNode):
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

    def get_post(self, user):
        return None

    def action(self, post, user, form):
        title = form.cleaned_data.get('title', '')
        type = form.cleaned_data.get('type', '')
        tags = form.cleaned_data.get('tags', '')
        content = form.cleaned_data.get('content', '')

        obj = Post.objects.create(content=content, title=title,
                                  author=user, type=type)

        # Needs to save root explicitly!
        obj.root.save()

        return obj
