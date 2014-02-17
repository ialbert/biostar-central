"""
Moderator views
"""
from biostar.apps.posts.models import Post
from django.conf import settings
from django.views.generic import  FormView
from biostar.apps.posts.auth import post_permissions
from django.shortcuts import render
from django.contrib import messages
from biostar import const
from braces.views import LoginRequiredMixin
from django import forms
from django.core.urlresolvers import reverse
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Fieldset, Submit, ButtonHolder
from django.http import HttpResponseRedirect

OPEN, CLOSE_OFFTOPIC, CLOSE_SPAM, DELETE, DUPLICATE = map(str, range(5))

class PostModForm(forms.Form):
    CHOICES = [
        (OPEN, "Open a closed or deleted post"),
        (DUPLICATE, "Close post, duplicated"),
        (CLOSE_OFFTOPIC, "Close post, off topic"),
        (CLOSE_SPAM, "Close post, inappropriate"),
        (DELETE, "Delete post"),
    ]

    action = forms.ChoiceField(choices=CHOICES, widget=forms.RadioSelect(), label="Select Action")

    comment = forms.CharField(required=False, max_length=200,
                             help_text="Enter a reason (required when closing or deleting, optional for duplicates)")

    dupe_title = forms.CharField(required=False, max_length=200,
                             help_text="Enter duplicated post number or title (required for duplicates, will link the duplicate)",
                             label="Duplicate number/title")

    def __init__(self, *args, **kwargs):
        pk = kwargs['pk']
        kwargs.pop('pk')
        super(PostModForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.error_text_inline = False
        self.helper.help_text_inline = True
        self.helper.form_action = reverse("post-moderation", kwargs=dict(pk=pk))

        self.helper.layout = Layout(
            Fieldset(
                'Select moderation option',
                'action',
                'comment',
                'dupe_title',
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )

class PostModeration(LoginRequiredMixin, FormView):
    model = Post
    template_name = "moderator-panel.html"
    context_object_name = "post"
    form_class = PostModForm

    def get_obj(self):
        pk = self.kwargs['pk']
        obj = Post.objects.get(pk=pk)
        return obj

    def get(self, request, *args, **kwargs):
        post = self.get_obj()
        post = post_permissions(request, post)
        if not post.is_editable:
            messages.warning(request, "You may not moderate this post")
            return HttpResponseRedirect(post.root.get_absolute_url())
        form = self.form_class(pk=post.id)
        context = dict(form=form, post=post)
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        user = request.user

        post = self.get_obj()
        post = post_permissions(request, post)

        if not post.is_editable:
            messages.warning(request, "You may not moderate this post")
            return HttpResponseRedirect(post.root.get_absolute_url())

        form = self.form_class(request.POST, pk=post.id)

        # Using a new query to update. Bypasses signals.
        query = Post.objects.filter(pk=post.id)
        if form.is_valid():
            get = form.cleaned_data.get
            action = get('action')
            if action == OPEN and not user.is_moderator:
                messages.error(request, "Only a moderator may open this post")
            elif action == OPEN:
                query.update(status=Post.OPEN)
                messages.success(request, "Opened %s" % post.title)
            elif action in (CLOSE_OFFTOPIC, CLOSE_SPAM, DUPLICATE):
                query.update(status=Post.CLOSED)
                messages.success(request, "Closed %s" % post.title)

                # Insert the comment or duplicated text
                comment = Post(content="Closed for some reason", type=Post.COMMENT, parent=post, author=user)
                comment.save()

            elif action == DELETE:
                children = Post.objects.filter(parent_id=post.id).exclude(pk=post.id)
                if children:
                    query.update(status=Post.DELETED)
                    messages.success(request, "Deleted %s" % post.title)
                else:
                    post.delete()
                    messages.success(request, "Destroyed %s" % post.title)
                    HttpResponseRedirect("/")
            else:
                messages.error(request, "Invalid action %s" % action)

        return HttpResponseRedirect(post.root.get_absolute_url())

    #    pid = request.GET.get("post_id")
    #    context = {}
    #    return render(request, self.template_name, context)