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


def post_title_validator(text):
    "Validates form input for tags"
    text = text.strip()
    if not text:
        raise ValidationError('Please enter a title')

    if len(text) < 10:
        raise ValidationError('The title is too short')

    words = text.split(" ")
    if len(words) < 3:
        raise ValidationError('Please have more than two words in the title.')


def post_id_validator(text):
    try:
        value = int(text)
    except ValueError, exc:
        raise ValidationError("The post_id must be an integer")

    try:
        post = Post.objects.get(pk=value)
    except ObjectDoesNotExist, exc:
        raise ValidationError("The post does not exist. Perhaps it was deleted")

        # Need to validate access to the post


class ContentForm(forms.Form):
    """
    Edit or create form for content: answers, comments
    """
    #action = forms.CharField(widget=forms.HiddenInput, required=True)
    #post_id = forms.IntegerField(validators=[post_id_validator], widget=forms.HiddenInput, required=True)
    content = forms.CharField(widget=forms.Textarea,
                              min_length=settings.MIN_POST_SIZE, max_length=settings.MAX_POST_SIZE,
                              initial="", required=True)

    #def __init__(self, *args, **kwargs):
    #    # Needs to know about the request to validate access.
    #    self.request = kwargs.pop('request', None)
    #    super(ContentForm, self).__init__(*args, **kwargs)

    '''
    def clean(self):
        cleaned_data = super(ContentForm, self).clean()
        post_id = cleaned_data.get('post_id', 0)

        post = Post.objects.filter(pk=post_id).select_related("group", "group__groupinfo").first()

        # The post may not exist anymore.
        if not post:
            raise ValidationError("Post does not exist. Perhaps it has been deleted.")

        # User must be in the group that created the thread this post belongs to.
        if not auth.read_access_post(user=self.request.user, post=post):
            raise ValidationError("Post may not be accessed by this user!")
    '''

def redirect(name):
    return HttpResponseRedirect(reverse(name))

class BaseNode(LoginRequiredMixin, FormView):
    form_class = ContentForm
    template_name = "edit_content.html"

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

    def action(self, post, user, content):
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
            self.post = self.action(post=post, user=user, content=content)
        except Exception, exc:
            logger.error(exc)
            messages.error(self.request, "Unable to create the post!")
            return redirect("home")

        return super(BaseNode, self).form_valid(form)

    def get_success_url(self):
        return self.post.get_absolute_url()


class NewNode(BaseNode):
    def action(self, post, user, content):
        post = Post.objects.create(parent=post, content=content, author=user)
        return post

class EditNode(BaseNode):

    def get_initial(self):
        try:
           post = self.get_post(user=self.request.user)
        except auth.AccessDenied, exc:
            logger.error(exc)
            return redirect("home")
        initial = dict(content=post.content)
        return initial

    def action(self, post, user, content):
        # Update the post
        post.content = content
        post.lastedit_user = user
        post.lastedit_date = auth.now()
        post.save()
        return post

