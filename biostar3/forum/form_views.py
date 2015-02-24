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

import logging

logger = logging.getLogger('biostar')


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


def parent_id_validator(text):
    try:
        value = int(text)
    except ValueError, exc:
        raise ValidationError("The parent_id must be an integer")
    try:
        parent = Post.objects.get(pk=value)
    except ObjectDoesNotExist, exc:
        raise ValidationError("The parent does not exist. Perhaps it was deleted")


def post_id_validator(text):
    try:
        value = int(text)
    except ValueError, exc:
        raise ValidationError("The post_id must be an integer")

    if value == 0:
        # This indicates a new entry into the database. Allow to proceed.
        return

    try:
        post = Post.objects.get(pk=value)
    except ObjectDoesNotExist, exc:
        raise ValidationError("The post does not exist. Perhaps it was deleted")

    # Need to validate access to the post


class NewContentForm(forms.Form):
    """
    For posts that have only content: answers, comments
    """
    parent_id = forms.IntegerField(validators=[parent_id_validator], widget=forms.HiddenInput, required=True)
    content = forms.CharField(widget=forms.Textarea,
                              min_length=settings.MIN_POST_SIZE, max_length=settings.MAX_POST_SIZE,
                              label="Enter your answer", initial="")


class NewContent(LoginRequiredMixin, FormView):
    form_class = NewContentForm
    template_name = "new_content.html"

    def form_valid(self, form):
        # Create the new post based on access rights.
        user = self.request.user
        parent_id = form.cleaned_data['parent_id']
        content = form.cleaned_data['content']
        try:
            parent = Post.objects.get(pk=parent_id)
            self.post = Post.objects.create(parent=parent, content=content, author=user,
                                            type=Post.ANSWER, lastedit_user=user)
        except KeyError, exc:
            logger.error(exc)
            self.post = None
        return super(NewContent, self).form_valid(form)

    def get_success_url(self):
        # Redirect to new post
        if self.post:
            print self.post.id
            return self.post.get_absolute_url()
        else:
            messages.error(self.request, "Unable to create the post")
            return reverse("home")

class EditContent(FormView):
    """
    Edits existing content.
    """
    form_class = NewContentForm
    template_name = "new_content.html"

    def get(self, request, *args, **kwargs):
        form_class = self.get_form_class()
        form = self.get_form(form_class)
        context = self.get_context_data(form=form)

        return self.render_to_response(context)


    def form_valid(self, form):
        self.obj = 123
        return super(EditContent, self).form_valid(form)

    def get_success_url(self):
        messages.info(self.request, "SUCCESS %s" % self.obj)
        return '/'