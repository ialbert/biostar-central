"""
Moderator views
"""
from biostar.apps.posts.models import Post
from biostar.apps.util import html

from django.conf import settings
from django.views.generic import FormView
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
        (DUPLICATE, "Duplicated, close"),
        (CLOSE_OFFTOPIC, "Closing, off topic"),
        (DELETE, "Delete post"),
    ]

    action = forms.ChoiceField(choices=CHOICES, widget=forms.RadioSelect(), label="Select Action")

    comment = forms.CharField(required=False, max_length=200,
                              help_text="Enter a reason (required when closing)")

    dupe = forms.CharField(required=False, max_length=200,
                           help_text="One or more duplicated post numbers, space or comma separated (required for duplicate closing).",
                           label="Duplicate number(s)")

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
                'dupe',
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )

    def clean(self):
        cleaned_data = super(PostModForm, self).clean()
        action = cleaned_data.get("action")
        comment = cleaned_data.get("comment")
        dupe = cleaned_data.get("dupe")

        if action == CLOSE_OFFTOPIC and not comment:
            raise forms.ValidationError("Unable to close. Please add a comment!")

        if action == DUPLICATE and not dupe:
            raise forms.ValidationError("Unable to close duplicate. Please fill in the post numbers")

        if dupe:
            dupe = dupe.replace(",", " ")
            dupes = dupe.split()[:5]
            cleaned_data['dupe'] = dupes

        return cleaned_data

class PostModeration(LoginRequiredMixin, FormView):
    model = Post
    template_name = "moderator_form.html"
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

        # Bail out on errors.
        if not form.is_valid():
            messages.error(request, "%s" % form.errors)
            return HttpResponseRedirect(post.root.get_absolute_url())

        get = form.cleaned_data.get

        # Using a new query to update. Bypasses signals.
        query = Post.objects.filter(pk=post.id)


        action = get('action')
        if action == OPEN and not user.is_moderator:
            messages.error(request, "Only a moderator may open this post")

        elif action == OPEN:
            query.update(status=Post.OPEN)
            messages.success(request, "Opened %s" % post.title)

        elif action in CLOSE_OFFTOPIC:
            query.update(status=Post.CLOSED)
            messages.success(request, "Closed %s" % post.title)
            content = html.render(name="messages/offtopic_posts.html", user=post.author, comment=get("comment"), post=post)
            comment = Post(content=content, type=Post.COMMENT, parent=post, author=user)
            comment.save()

        elif action == DUPLICATE:
            query.update(status=Post.CLOSED)
            posts = Post.objects.filter(id__in=get("dupe"))
            content = html.render(name="messages/duplicate_posts.html", user=post.author, comment=get("comment"), posts=posts)
            comment = Post(content=content, type=Post.COMMENT, parent=post, author=user)
            comment.save()

        elif action == DELETE:

            children = Post.objects.filter(parent_id=post.id).exclude(pk=post.id)

            #
            # Delete vs destroy
            #
            # Posts with children or older than some value can only be deleted not destroyed
            delete_only = children or post.age_in_days > 7 or post.vote_count > 1

            if delete_only:
                query.update(status=Post.DELETED)
                messages.success(request, "Deleted %s" % post.title)

            else:
                post.delete()
                messages.success(request, "Destroyed %s" % post.title)
                return HttpResponseRedirect("/")

        else:
            messages.error(request, "Invalid action %s" % action)

        return HttpResponseRedirect(post.root.get_absolute_url())

    #    pid = request.GET.get("post_id")
    #    context = {}
    #    return render(request, self.template_name, context)