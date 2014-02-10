# Create your views here.
from django.shortcuts import render_to_response
from django.views.generic import TemplateView, DetailView, ListView, FormView, UpdateView
from .models import Post
from django import forms
from django.core.urlresolvers import reverse
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Fieldset, Submit, ButtonHolder
from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.contrib import messages
from . import auth
from braces.views import LoginRequiredMixin
from datetime import datetime
from django.utils.timezone import utc

class LongForm(forms.Form):
    FIELDS = "title content".split()

    title = forms.CharField()
    content = forms.CharField(widget=forms.Textarea, required=False)

    def __init__(self, *args, **kwargs):
        super(LongForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.layout = Layout(
            Fieldset(
                'Post',
                'title',
                'content',
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )


class ShortForm(forms.Form):
    FIELDS = "content"

    title = forms.CharField()
    content = forms.CharField(widget=forms.Textarea, required=False)

    def __init__(self, *args, **kwargs):
        super(ShortForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.layout = Layout(
            Fieldset(
                'Post',
                'content',
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )


class NewPost(LoginRequiredMixin, FormView):
    """
    Creates a new post.
    """
    form_class = LongForm
    template_name = "post-edit.html"

    def get(self, request, *args, **kwargs):
        initial = {}

        # The parent id.
        pid = int(self.kwargs['pid'])
        form_class = ShortForm if pid else LongForm
        form = form_class(initial=initial)
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):
        user = request.user
        pid = int(self.kwargs['pid'])

        # Posts with a parent are not toplevel
        form_class = ShortForm if pid else LongForm

        # Validating the form.
        form = form_class(request.POST)
        if not form.is_valid():
            return render(request, self.template_name, {'form': form})

        # Valid forms start here.
        data = form.cleaned_data

        # Set the parent if it exists
        parent = Post.objects.get(pk=pid) if pid else None

        # Set the title based on the level of the object.
        if parent:
            title = parent.title
        else:
            title = data['title']

        # Create a new post.
        post = Post(
            title=title, content=data['content'], author=user, type=Post.FORUM,
            parent=parent,
        )

        messages.success(request, "Post created")
        post.save()

        return HttpResponseRedirect(post.get_absolute_url())


class EditPost(LoginRequiredMixin, FormView):
    """
    Edits an existing post.
    """

    # The template_name attribute must be specified in the calling apps.
    template_name = "post-edit.html"
    form_class = LongForm

    def get(self, request, *args, **kwargs):
        initial = {}

        pk = int(self.kwargs['pk'])
        post = Post.objects.get(pk=pk)
        post = auth.post_permissions(request=request, post=post)

        # Check and exit if not a valid edit.
        if not post.is_editable:
            messages.error(request, "This user may not modify the post")
            return HttpResponseRedirect(reverse("home"))

        initial = dict(title=post.title, content=post.content)

        form_class = LongForm if post.is_toplevel else ShortForm
        form = form_class(initial=initial)
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):

        pk = int(self.kwargs['pk'])
        post = Post.objects.get(pk=pk)
        post = auth.post_permissions(request=request, post=post)

        # Check and exit if not a valid edit.
        if not post.is_editable:
            messages.error(request, "This user may not modify the post")
            return HttpResponseRedirect(reverse("home"))

        # Posts with a parent are not toplevel
        form_class = LongForm if post.is_toplevel else ShortForm

        form = form_class(request.POST)
        if not form.is_valid():
            # Invalid form submission.
            return render(request, self.template_name, {'form': form})

        # Valid forms start here.
        data = form.cleaned_data

        # Set the form attributes.
        for field in form_class.FIELDS:
            setattr(post, field, data[field])

        # Update the last editing user.
        post.lastedit_user = request.user
        post.lastedit_date = datetime.utcnow().replace(tzinfo=utc)
        post.save()
        messages.success(request, "Post updated")

        return HttpResponseRedirect(post.get_absolute_url())

    def get_success_url(self):
        return reverse("user-details", kwargs=dict(pk=self.kwargs['pk']))

