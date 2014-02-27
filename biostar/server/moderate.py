"""
Moderator views
"""
from biostar.apps.posts.models import Post
from biostar.apps.posts.auth import post_permissions

from biostar.apps.users.models import User
from biostar.apps.users.auth import user_permissions

from biostar.apps.util import html

from django.conf import settings
from django.views.generic import FormView

from django.shortcuts import render
from django.contrib import messages
from biostar import const
from braces.views import LoginRequiredMixin
from django import forms
from django.core.urlresolvers import reverse
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Fieldset, Submit, ButtonHolder
from django.http import HttpResponseRedirect
from django.db.models import Q, F

OPEN, CLOSE_OFFTOPIC, CLOSE_SPAM, DELETE, DUPLICATE, MOVE_TO_COMMENT, MOVE_TO_ANSWER = map(str, range(7))


class PostModForm(forms.Form):
    CHOICES = [
        (OPEN, "Open a closed or deleted post"),
        (MOVE_TO_ANSWER, "Move post to an answer"),
        (MOVE_TO_COMMENT, "Move post to a comment on the top level post"),
        (DUPLICATE, "Duplicated, close"),
        (CLOSE_OFFTOPIC, "Closing"),
        (DELETE, "Delete post"),
    ]

    action = forms.ChoiceField(choices=CHOICES, widget=forms.RadioSelect(), label="Select Action")

    comment = forms.CharField(required=False, max_length=200,
                              help_text="Enter a reason (required when closing). This will be inserted into a template comment.")

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
    template_name = "post_moderation_form.html"
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

        # The default return url
        response = HttpResponseRedirect(post.root.get_absolute_url())

        if not post.is_editable:
            messages.warning(request, "You may not moderate this post")
            return response

        # Initialize the form class.
        form = self.form_class(request.POST, pk=post.id)

        # Bail out on errors.
        if not form.is_valid():
            messages.error(request, "%s" % form.errors)
            return response

        # A shortcut to the clean form data.
        get = form.cleaned_data.get

        # These will be used in updates, will bypasses signals.
        query = Post.objects.filter(pk=post.id)
        root  = Post.objects.filter(pk=post.root_id)

        action = get('action')
        if action == OPEN and not user.is_moderator:
            messages.error(request, "Only a moderator may open a post")
            return response

        if action == MOVE_TO_ANSWER and post.type == Post.COMMENT:
            # This is a valid action only for comments.
            messages.success(request, "Moved post to answer")
            query.update(type=Post.ANSWER, parent=post.root)
            root.update(reply_count=F("reply_count") + 1)
            return response

        if action == MOVE_TO_COMMENT and post.type == Post.ANSWER:
            # This is a valid action only for answers.
            messages.success(request, "Moved post to answer")
            query.update(type=Post.COMMENT, parent=post.root)
            root.update(reply_count=F("reply_count") - 1)
            return response

        # Some actions are valid on top level posts only.
        if action in (CLOSE_OFFTOPIC, DUPLICATE, OPEN) and not post.is_toplevel:
            messages.warning(request, "You can only close or open a top level post")
            return response

        if action == OPEN:
            query.update(status=Post.OPEN)
            messages.success(request, "Opened post: %s" % post.title)
            return response

        if action in CLOSE_OFFTOPIC:
            query.update(status=Post.CLOSED)
            messages.success(request, "Closed post: %s" % post.title)
            content = html.render(name="messages/offtopic_posts.html", user=post.author, comment=get("comment"), post=post)
            comment = Post(content=content, type=Post.COMMENT, parent=post, author=user)
            comment.save()
            return response

        if action == DUPLICATE:
            query.update(status=Post.CLOSED)
            posts = Post.objects.filter(id__in=get("dupe"))
            content = html.render(name="messages/duplicate_posts.html", user=post.author, comment=get("comment"), posts=posts)
            comment = Post(content=content, type=Post.COMMENT, parent=post, author=user)
            comment.save()
            return response

        if action == DELETE:

            # Delete means to mark a post deleted. Remove means to delete the post from the database.
            # Posts with children or older than some value can only be deleted not removed
            # Banning a user will remove their posts.

            children = Post.objects.filter(parent_id=post.id).exclude(pk=post.id)
            delete_only = children or post.age_in_days > 7 or post.vote_count > 1 or (post.author != user)

            if delete_only:
                # Deleted posts can be undeleted by re-opening them.
                query.update(status=Post.DELETED)
                messages.success(request, "Deleted post: %s" % post.title)
                return response
            else:
                # This will remove the post. Redirect depends on the level of the post.
                url = "/" if post.is_toplevel else post.parent.get_absolute_url()
                post.delete()
                messages.success(request, "Removed post: %s" % post.title)
                return HttpResponseRedirect(url)

        # By this time all actions should have been performed
        messages.warning(request, "That seems to be an invalid action for that post. \
                It is probably ok! Actions may be shown even when not valid.")
        return response

class UserModForm(forms.Form):
    CHOICES = [
        (User.NEW_USER, "Reinstate as new user"),
        (User.TRUSTED, "Reinstate as trusted user"),
        (User.SUSPENDED, "Suspend user"),
        (User.BANNED, "Ban user"),
    ]

    action = forms.ChoiceField(choices=CHOICES, widget=forms.RadioSelect(), label="Select Action")

    def __init__(self, *args, **kwargs):
        pk = kwargs['pk']
        kwargs.pop('pk')
        super(UserModForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.error_text_inline = False
        self.helper.help_text_inline = True
        self.helper.form_action = reverse("user-moderation", kwargs=dict(pk=pk))

        self.helper.layout = Layout(
            Fieldset(
                'Select action',
                'action',
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )


class UserModeration(LoginRequiredMixin, FormView):
    model = Post
    template_name = "user_moderation_form.html"
    context_object_name = "user"
    form_class = UserModForm

    def get_obj(self):
        pk = self.kwargs['pk']
        obj = User.objects.get(pk=pk)
        return obj

    def get(self, request, *args, **kwargs):
        user = request.user
        target = self.get_obj()
        target = user_permissions(request, target)

        form = self.form_class(pk=target.id)
        context = dict(form=form, target=target)
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        user = request.user

        target = self.get_obj()
        target = user_permissions(request, target)

        # The response after the action
        response = HttpResponseRedirect(target.get_absolute_url())

        if not user.is_moderator or not target.is_editable:
            messages.warning(request, "Current user does not have sufficient moderator privileges")
            return response

        if target.is_administrator:
            messages.warning(request, "Cannot moderate an administrator")
            return response

        if user == target:
            messages.warning(request, "Cannot moderate yourself")
            return response

        form = self.form_class(request.POST, pk=target.id)
        if not form.is_valid():
            messages.error(request, "Invalid user modification action")
            return response

        action = int(form.cleaned_data['action'])

        if action == User.BANNED and not user.is_administrator:
            messages.info(request, "Only administrators can ban users")
            return response

        if action == User.BANNED and user.is_administrator:
            query = Post.objects.filter(author=target, type__in=Post.TOP_LEVEL, vote_count=0)
            count = query.count()
            query.delete()
            messages.success(request, "User banned, %s posts removed" % count)
            return response

        # Apply the new status
        User.objects.filter(pk=target.id).update(status=action)

        messages.success(request, 'Moderation completed')
        return response

