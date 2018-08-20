from django import forms
from .models import Post
from django.core.exceptions import ValidationError
from django.db.models import F
from django.conf import settings

from biostar.engine.models import Project
from . import  models, auth
from pagedown.widgets import PagedownWidget

# Share logger with models
logger = models.logger


def english_only(text):


    try:
        text.encode('ascii')
    except Exception:
        raise ValidationError('Title may only contain plain text (ASCII) characters')


def valid_title(text):
    "Validates form input for tags"
    text = text.strip()
    if not text:
        raise ValidationError('Please enter a title')

    words = text.split()
    if len(words) < 3:
        raise ValidationError('More than two words please.')


def valid_tag(text):
    "Validates form input for tags"

    words = text.split(",")
    if len(words) > 5:
        raise ValidationError('You have too many tags (5 allowed)')


class PostLongForm(forms.Form):

    def __init__(self, project=None, filter_func=lambda x: x, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.project = project

        POST_CHOICES = [(Post.QUESTION, "Question"),
                        (Post.JOB, "Job Ad"),
                        (Post.TUTORIAL, "Tutorial"), (Post.TOOL, "Tool"),
                        (Post.FORUM, "Forum"), (Post.NEWS, "News"),
                        (Post.BLOG, "Blog"), (Post.PAGE, "Page")]

        # Pass a filtering function to customize choices between sites.
        choices = list(filter(filter_func, POST_CHOICES))
        self.fields["post_type"] = forms.ChoiceField(label="Post Type", choices=choices,
                                                     help_text="Select a post type: Question, Forum, Job, Blog")

    title = forms.CharField(
        label="Post Title",
        max_length=200, min_length=2, validators=[valid_title, english_only],
        help_text="Descriptive titles promote better answers.")

    tag_val = forms.CharField(
        label="Post Tags", max_length=50,
        required=False, validators=[valid_tag],
        help_text="Choose one or more tags to match the topic. To create a new tag just type it in comma seperated.",
    )

    content = forms.CharField(widget=PagedownWidget(template="widgets/pagedown.html"), validators=[english_only],
                              min_length=10, max_length=15000,
                              label="Enter your post below")

    def save(self, author=None):
        data = self.cleaned_data.get

        title = data('title')
        content = data('content')
        post_type = int(data('post_type'))
        tag_val = data('tag_val')

        post = auth.create_post(title=title, content=content, post_type=post_type,
                                tag_val=tag_val, author=author, project=self.project)

        return post


class SubsForm(forms.Form):

    choices = models.Subscription.MESSAGING_CHOICES

    subtype = forms.IntegerField(widget=forms.Select(choices=choices))

    def __init__(self, user, post, *args, **kwargs):

        self.user = user
        self.post = post

        super(SubsForm, self).__init__(*args, **kwargs)

    def save(self):

        sub_type = self.cleaned_data["subtype"]

        sub = auth.create_sub(post=self.post, user=self.user, sub_type=sub_type)

        return sub


class PostShortForm(forms.Form):

    content = forms.CharField(widget=PagedownWidget(template="widgets/pagedown.html"),
                              min_length=2, max_length=5000)

    parent_uid = forms.CharField(widget=forms.HiddenInput(), min_length=2, max_length=5000)
    project_uid = forms.CharField(widget=forms.HiddenInput(), min_length=2, max_length=5000,
                                  required=False)
    redir_url = forms.CharField(widget=forms.HiddenInput(), min_length=2, max_length=5000,
                                  required=True)

    def save(self, author, post_type=Post.ANSWER):
        data = self.cleaned_data

        parent = Post.objects.filter(uid=data.get("parent_uid")).first()
        project = Project.objects.filter(uid=data.get("project_uid")).first()
        auth.create_post(title=parent.title,
                          parent=parent,
                          author=author,
                          content=data.get("content"),
                          post_type=post_type,
                          project=project,
                          sub_to_root=True
                          )

        return data.get("redir_url", "/")

    #def clean(self):
    #    cleaned_data = super(PostShortForm, self).clean()
    #    return








