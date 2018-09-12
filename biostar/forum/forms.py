import mistune
from pagedown.widgets import PagedownWidget
from django import forms
from .models import Post
from django.core.exceptions import ValidationError
from django.db.models import F
from django.conf import settings
from biostar.engine.models import Project

from . import  models, auth, util

from .const import *
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

    def __init__(self, project=None, filter_func=lambda x: x, post=None, user=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.project = project
        self.post = post
        self.user = user

        POST_CHOICES = [(Post.QUESTION, "Question"),
                        (Post.JOB, "Job Ad"),
                        (Post.TUTORIAL, "Tutorial"), (Post.TOOL, "Tool"),
                        (Post.FORUM, "Forum"), (Post.NEWS, "News"),
                        (Post.BLOG, "Blog"), (Post.PAGE, "Page")]

        # Pass a filtering function to customize choices between sites.
        choices = list(filter(filter_func, POST_CHOICES))
        inital_title = inital_tags = inital_content = ""
        inital_type = Post.FORUM
        if self.post:
            inital_title = self.post.title
            inital_tags = self.post.tag_val
            inital_content = self.post.html
            inital_type = self.post.type

        self.fields["post_type"] = forms.ChoiceField(label="Post Type", choices=choices,
                                                     help_text="Select a post type: Question, Forum, Job, Blog",
                                                     initial=inital_type)
        self.fields["title"] = forms.CharField(label="Post Title", max_length=200, min_length=2, initial=inital_title,
                                               validators=[valid_title, english_only],
                                               help_text="Descriptive titles promote better answers.")
        self.fields["tag_val"] = forms.CharField(label="Post Tags", max_length=50, required=False, validators=[valid_tag],
                                                 help_text="""Choose one or more tags to match the topic. 
                                                        To create a new tag just type it in comma seperated.""",
                                                 initial=inital_tags)
        self.fields["content"] = forms.CharField(widget=PagedownWidget(template="widgets/pagedown.html"), validators=[english_only],
                                                 min_length=10, max_length=15000, initial=inital_content,
                                                 label="Enter your post below")

    def save(self, author=None, edit=False):
        data = self.cleaned_data.get

        title = data('title')
        html = data("content")
        content = util.strip_tags(html)
        post_type = int(data('post_type'))
        tag_val = data('tag_val')

        if edit:
            self.post.title = title
            self.post.content = content
            self.post.type = post_type
            self.post.html = html
            self.post.save()
            # Triggers another save
            self.post.add_tags(text=tag_val)
        else:
            self.post = auth.create_post(title=title, content=content, post_type=post_type,
                                         tag_val=tag_val, author=author, project=self.project)

        return self.post

    def clean(self):

        if self.post and self.user != self.post.author:
            if self.user.profile.is_manager or self.user.profile.is_moderator:
                pass
            else:
                raise forms.ValidationError("Only the author or a moderator can edit a post.")


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

    def __init__(self, user=None, post=None, *args, **kwargs):

        self.user = user
        self.post = post

        super().__init__(*args, **kwargs)
        inital_content = "" if not self.post else self.post.html
        self.fields["content"] = forms.CharField(widget=PagedownWidget(template="widgets/pagedown.html"),
                                                 min_length=2, max_length=5000,
                                                 initial=inital_content)

    parent_uid = forms.CharField(widget=forms.HiddenInput(), min_length=2, max_length=5000)
    project_uid = forms.CharField(widget=forms.HiddenInput(), min_length=2, max_length=5000,
                                  required=False)

    def save(self, author=None, post_type=Post.ANSWER, edit=False):
        data = self.cleaned_data
        html = data.get("content")
        project = data.get("project_uid")
        parent = data.get("parent_uid")
        content = util.strip_tags(html)

        if edit:
            self.post.html = html
            self.post.content = content
            self.post.save()
        else:
            parent = Post.objects.get_all(uid=parent).first()
            project = Project.objects.filter(uid=project).first()

            self.post = auth.create_post(title=parent.root.title,
                              parent=parent,
                              author=author,
                              content=content,
                              post_type=post_type,
                              project=project,
                              sub_to_root=True,
                              )
        return self.post

    #def clean(self):
    #    cleaned_data = super(PostShortForm, self).clean()
    #    return


class PostModForm(forms.Form):

    CHOICES = [
        (BUMP_POST, "Bump a post"),
        (MOD_OPEN, "Open a closed or deleted post"),
        (TOGGLE_ACCEPT, "Toggle accepted status"),
        (MOVE_TO_ANSWER, "Move post to an answer"),
        (MOVE_TO_COMMENT, "Move post to a comment on the top level post"),
        (DUPLICATE, "Duplicated post (top level)"),
        (CROSSPOST, "Cross posted at other site"),
        (CLOSE_OFFTOPIC, "Close post (top level)"),
        (DELETE, "Delete post"),
    ]

    action = forms.ChoiceField(choices=CHOICES, widget=forms.RadioSelect(), label="Select Action")

    comment = forms.CharField(required=False, max_length=200,
                              help_text="Enter a reason (required when closing, crosspost). This will be inserted into a template comment.")

    dupe = forms.CharField(required=False, max_length=200,
                           help_text="One or more duplicated post numbers, space or comma separated (required for duplicate closing).",
                           label="Duplicate number(s)")

    def __init__(self, post, request, user, *args, **kwargs):
        self.post = post
        self.user = user
        self.request = request
        super(PostModForm, self).__init__(*args, **kwargs)

    def save(self):

        cleaned_data = self.cleaned_data
        action = int(cleaned_data.get("action"))
        comment = cleaned_data.get("comment")
        dupe = cleaned_data.get("dupe")

        url = auth.moderate_post(post=self.post, request=self.request,
                                 action=action, comment=comment, dupes=dupe)
        return url

    def clean(self):
        cleaned_data = super(PostModForm, self).clean()
        action = int(cleaned_data.get("action"))
        comment = cleaned_data.get("comment")
        dupe = cleaned_data.get("dupe")

        if not (self.user.profile.is_moderator or self.user.profile.is_manager):
            raise forms.ValidationError( "Only a moderator/manager may perform these actions")

        if action in (CLOSE_OFFTOPIC, DUPLICATE, BUMP_POST) and not self.post.is_toplevel:
            raise forms.ValidationError("You can only perform these actions to a top-level post")
        elif action in (TOGGLE_ACCEPT, MOVE_TO_COMMENT) and self.post.type != Post.ANSWER:
            raise forms.ValidationError("You can only perform these actions to an answer.")
        elif action == MOVE_TO_ANSWER and self.post.type != Post.COMMENT:
            raise forms.ValidationError("You can only perform these actions to a comment.")
        if action == CLOSE_OFFTOPIC and not comment:
            raise forms.ValidationError("Unable to close. Please add a comment!")

        if action == CROSSPOST and not comment:
            raise forms.ValidationError("Please add URL into the comment!")

        if action == DUPLICATE and not dupe:
            raise forms.ValidationError("Unable to close duplicate. Please fill in the post numbers")

        if dupe:
            dupe = dupe.replace(",", " ")
            dupes = dupe.split()[:5]
            cleaned_data['dupe'] = dupes

        return cleaned_data




