from __future__ import absolute_import, division, print_function, unicode_literals

import logging
from datetime import timedelta

from django import forms
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.core.exceptions import ValidationError
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect
from django.db.models import Q
from django.db import transaction
from django.contrib.auth import get_user_model
from .models import Post, Profile, right_now, Site
from . import auth, cache

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

    file = forms.FileField(widget=forms.ClearableFileInput, required=False, label="File")


def get_post_form(request):
    sites = Site.objects.all().order_by("id")
    site_choices = [(site.id, site.name) for site in sites]

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
            (Post.QUESTION, "Question"),
            (Post.NEWS, "News"),
            (Post.FORUM, "Forum"),
            (Post.JOB, "Job Ad"),
            (Post.PAGE, "Page"),
        ])

        site = forms.TypedChoiceField(coerce=int, choices=site_choices, initial=request.site.id)

    return PostForm


@transaction.atomic
def post_create(request, parent=None, post_type=None, action='', top_level=False):
    """
    This view creates nodes. Is not called directly from the web only through
    other functions that prefill parameters.
    """
    user = request.user
    template_name = "post_edit.html"

    # Pick the right form to operate on.
    if top_level:
        form_class = get_post_form(request)
    else:
        form_class = ContentForm

    if request.method == "GET":
        # This will render the initial form for the user.
        form = form_class
        context = dict(form=form, action=action)
        return render(request, template_name, context)

    if request.method == "POST":
        # Data is being submitted
        form = form_class(request.POST)

        # Attempt to detect duplicated submissions
        content = request.POST.get("content")

        recently = right_now() - timedelta(minutes=10)
        if Post.objects.filter(content=content, creation_date__gt=recently):
            form.add_error("content", "Duplicated submission? There is a recent post with identical content!")

        if post_type is not None and parent is None:
            form.add_error("type", "Top level post may not have a parent.")

        # Handle the file upload
        file = request.FILES.get('file')

        if file and file.size > settings.MAX_UPLOAD_SIZE:
            maxsize = settings.MAX_UPLOAD_SIZE / 1024.0 / 1024
            form.add_error("content", "The uploaded file is too large! only %4.1f MB allowed" % maxsize)

        if not form.is_valid():
            # Form data came but not valid.
            context = dict(form=form, action=action)
            return render(request, template_name, context)


        # The form is valid create the post based on the form.
        if post_type is None:
            post = auth.create_toplevel_post(user=user, data=form.cleaned_data, file=file)
        else:
            content = form.cleaned_data['content']
            post = auth.create_content_post(content=content, post_type=post_type, user=user, parent=parent, file=file)

        return redirect(post.get_absolute_url())

    messages.error(request, "Unsupported request type")
    return redirect("home")


@login_required
def create_toplevel_post(request):
    "A new toplevel post"
    action = reverse("new_post")
    return post_create(request=request, parent=None, post_type=None, action=action, top_level=True)


@login_required
@auth.content_create
def create_answer(request, pk, parent=None):
    action = reverse("new_answer", kwargs=dict(pk=parent.id))
    return post_create(request=request, parent=parent, post_type=Post.ANSWER, action=action, top_level=False)


@login_required
@auth.content_create
def create_comment(request, pk, parent=None):
    action = reverse("new_comment", kwargs=dict(pk=parent.id))
    return post_create(request=request, parent=parent, post_type=Post.COMMENT, action=action, top_level=False)


@login_required
@auth.post_edit
@transaction.atomic
def post_edit(request, pk, post=None, user=None):
    """
    This view updates posts.
    """

    template_name = "post_edit.html"
    action = reverse("post_edit", kwargs=dict(pk=post.id))

    if post.is_toplevel:
        # Different forms are chosen based on post type.
        # A form with title, type and tags.
        form_class = get_post_form()
        tags = ", ".join(post.tags.names())
        initial = dict(content=post.content, title=post.title, tags=tags, site=post.site.id, type=post.type,
                       file=post.file)
    else:
        # Content only: answers and comments.
        form_class = ContentForm
        initial = dict(content=post.content, file=post.file)

    if request.method == "GET":
        # Get methods get the form and return.
        form = form_class(initial=initial)
        context = dict(form=form, action=action)
        return render(request, template_name, context)

    if request.method == "POST":
        # This is a form submission with incoming parameters.
        form = form_class(request.POST)

        # Handle the file upload allowing file removal.
        file = request.FILES.get('file')

        if file and file.size > settings.MAX_UPLOAD_SIZE:
            maxsize = settings.MAX_UPLOAD_SIZE / 1024.0 / 1024
            form.add_error("content", "The uploaded file is too large! only %4.1f MB allowed" % maxsize)

        if not form.is_valid():
            # Invalid form, bail out with error messaged.
            context = dict(form=form, action=action)
            return render(request, template_name, context)

        # The data is valid update the post and return the view to it.
        get = lambda word: form.cleaned_data.get(word, '')

        post.content = get('content')
        post.lastedit_user = user
        post.lastedit_date = right_now()

        clear = request.POST.get("file-clear")

        if file or clear:
            post.file = file

        # Extra information to be saved for toplevel posts.
        if post.is_toplevel:
            post.site_id = get('site')
            post.title = get('title')
            post.type = get('type')
            tags = get('tags')
            tags = auth.tag_split(tags)
            # Must explicitly set the new tags.
            post.tags.set(*tags)

        post.save()

    return redirect(post.get_absolute_url())


# Full lenght text input widget.
text_input = lambda: forms.TextInput(attrs={"class": "u-full-width"})


class UserProfileForm(forms.Form):
    name = forms.CharField(max_length=100, widget=text_input())
    handle = forms.CharField(max_length=100, widget=text_input())
    email = forms.CharField(max_length=150, widget=text_input())
    website = forms.CharField(max_length=150, required=False, widget=text_input())
    twitter_id = forms.CharField(max_length=150, label="Twitter ID", required=False, widget=text_input())
    scholar = forms.CharField(max_length=150, label="Google Scholar ID", required=False, widget=text_input())
    location = forms.CharField(max_length=150, required=False, widget=text_input())
    info = forms.CharField(widget=forms.Textarea, required=False, max_length=3000)
    shortcuts_text = forms.CharField(widget=forms.Textarea, required=False, max_length=500)

    def clean_email(self):
        """
        Can only set email to one that does not already exists for a different user
        """
        email = self.cleaned_data['email']
        user = User.objects.filter(email=email).first()
        if user and user != self.target:
            raise forms.ValidationError("This email already exists in the system for a different user!\
                    Contact us if you need to merge accounts.")
        return email

    def __init__(self, target, *args, **kwargs):
        self.target = target
        super(UserProfileForm, self).__init__(*args, **kwargs)


@login_required
@auth.valid_user
def user_edit(request, pk, target=None):
    template_name = "user_edit.html"

    if request.method not in ("GET", "POST"):
        messages.error(request, "Invalid request method")
        return redirect(target.get_absolute_url())

    if request.method == "GET":
        initial = dict(
            name=target.name,
            handle=target.handle,
            email=target.email,
            location=target.profile.location,
            website=target.profile.website,
            twitter_id=target.profile.twitter_id,
            scholar=target.profile.scholar,
            info=target.profile.info,
            shortcuts_text=target.profile.shortcuts_text,

        )
        form = UserProfileForm(target, initial=initial)
        context = dict(form=form, target=target)
        return render(request, template_name, context)

    form = UserProfileForm(target, data=request.POST)

    if not form.is_valid():
        context = dict(form=form, target=target)
        return render(request, template_name, context)

    # The form is valid. Process it.
    get = lambda key: form.cleaned_data.get(key, '')

    # Need to update both the user and the profile.
    User.objects.filter(pk=target.id).update(
        name=get("name"), email=get("email"), handle=get('handle'),
    )

    profile = Profile.objects.filter(user__id=target.id).first()
    for field in "info website scholar location twitter_id shortcuts_text scholar".split():
        setattr(profile, field, get(field))

    # Trigger save.
    profile.save()

    return redirect(target.get_absolute_url())


def suggest_feed(request):
    """
    Small form to suggest an RSS feed.
    """
    template_name = "planet_addfeed.html"

    add_feed = redirect("planet_addfeed")

    if request.method == "GET":
        context = dict()
        return render(request, template_name, context)
    else:
        feed = request.POST.get('feed', '').strip()
        recaptcha_response = request.POST.get('g-recaptcha-response')

        if not auth.valid_captcha(request):
            # Captcha not valid.
            return add_feed

    messages.info("Suggestion has been recorded. An admin will review the feed.")
    context = dict()
    return add_feed