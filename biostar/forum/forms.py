import mistune
from difflib import SequenceMatcher
from pagedown.widgets import PagedownWidget
from django import forms
from .models import Post
from django.core.exceptions import ValidationError
from django.db.models import F
from django.conf import settings
from biostar.engine.models import Project
from biostar.accounts.models import User
from biostar.forum.awards import *
from biostar.message.tasks import send_message, send_subs_msg, parse_mentioned_users, parse_mention_msg
from biostar.message.models import Message
from biostar.forum import models, auth, util

from .const import *
# Share logger with models
logger = models.logger


def send(old_content, new_content, post):
    """
    See if there is any change and send
    notifications and subscriptions messages
    """
    # Get rid of white spaces and tabs by splitting into lists.
    source_list = old_content.strip().split()
    new_list = new_content.strip().split()
    diffobj = SequenceMatcher(a=source_list, b=new_list)

    # There is a change detected
    change = diffobj.ratio() != 1

    # Do nothing when no change is detected.
    if not change:
        return

    # Get the sender for messages from biostar
    sender = User.objects.filter(is_superuser=True).first()

    # Send message to subscribed users.
    send_subs_msg(post=post)

    # Get mentioned users from new content
    new_mentioned_users = parse_mentioned_users(content=new_content)
    # Get old mentioned users and exclude them from these round of messages.
    old_mentioned_users = parse_mentioned_users(content=old_content).values("id")

    # Exclude old mentioned users from new ones.
    ment_users = new_mentioned_users.exclude(id__in=old_mentioned_users)

    # Parse the mentioned message
    ment_body, ment_subject, _ = parse_mention_msg(post=post)

    # Send the mentioned notifications
    send_message(source=Message.MENTIONED, subject=ment_subject, body=ment_body,
                 rec_list=ment_users, sender=sender)
    return


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
                        (Post.BLOG, "Blog")]

        # Pass a filtering function to customize choices between sites.
        choices = list(filter(filter_func, POST_CHOICES))
        inital_title = inital_tags = inital_content = ""
        inital_type = Post.FORUM
        if self.post:
            inital_title = self.post.title
            inital_tags = self.post.tag_val
            inital_content = self.post.content
            inital_type = self.post.type

        self.fields["post_type"] = forms.IntegerField(label="Post Type", widget=forms.Select(choices=choices),
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
        data = self.cleaned_data

        title = data.get('title')
        content = data.get("content")
        html = auth.parse_html(content)
        post_type = data.get('post_type')
        tag_val = data.get('tag_val')

        if edit:

            self.post.title = title
            old_content = self.post.content
            self.post.content = content
            self.post.type = post_type
            self.post.html = html
            self.post.tag_val = tag_val
            self.post.save()
            send(old_content=old_content, new_content=self.post.content, post=self.post)
            # Triggers another save
            #self.post.add_tags(text=tag_val)
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

        sub = models.Subscription.objects.filter(post=self.post, user=self.user).first()

        if sub and (sub.type == sub_type):
            return sub

        sub = auth.create_sub(post=self.post, user=self.user, sub_type=sub_type)

        return sub


class PostShortForm(forms.Form):

    def __init__(self, user=None, post=None, *args, **kwargs):

        self.user = user
        self.post = post

        super().__init__(*args, **kwargs)
        inital_content = "" if not self.post else self.post.content
        self.fields["content"] = forms.CharField(widget=PagedownWidget(template="widgets/pagedown.html"),
                                                 min_length=2, max_length=5000,
                                                 initial=inital_content)

    parent_uid = forms.CharField(widget=forms.HiddenInput(), min_length=2, max_length=5000)
    project_uid = forms.CharField(widget=forms.HiddenInput(), min_length=2, max_length=5000,
                                  required=False)

    def save(self, author=None, post_type=Post.ANSWER, edit=False):
        data = self.cleaned_data
        project = data.get("project_uid")
        parent = data.get("parent_uid")
        content = data.get("content")
        html = auth.parse_html(content)

        if edit:
            self.post.html = html
            self.post.content = content
            self.post.save()
        else:
            parent = Post.objects.filter(uid=parent).first()
            project = Project.objects.filter(uid=project).first()

            self.post = auth.create_post(title=parent.root.title,
                                         parent=parent,
                                         author=author,
                                         content=content,
                                         post_type=post_type,
                                         project=project,
                                         sub_to_root=True)
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
        (DELETE, "Delete post"),
    ]

    action = forms.IntegerField(widget=forms.RadioSelect(choices=CHOICES), label="Select Action")

    dupe = forms.CharField(required=False, max_length=200,
                           help_text="""One or more duplicated post numbers, 
                                       comma separated (required for duplicate closing).
                                       """,
                           label="Duplicate number(s)")

    def __init__(self, post, request, user, *args, **kwargs):
        self.post = post
        self.user = user
        self.request = request
        super(PostModForm, self).__init__(*args, **kwargs)

    def save(self):

        cleaned_data = self.cleaned_data
        action = cleaned_data.get("action")
        dupe = cleaned_data.get("dupe")

        url = auth.moderate_post(post=self.post, request=self.request,
                                 action=action, dupes=dupe)
        return url

    def clean(self):
        cleaned_data = super(PostModForm, self).clean()
        action = cleaned_data.get("action")
        dupe = cleaned_data.get("dupe")

        if not (self.user.profile.is_moderator or self.user.profile.is_manager):
            raise forms.ValidationError("Only a moderator/manager may perform these actions")

        if action in (DUPLICATE, BUMP_POST) and not self.post.is_toplevel:
            raise forms.ValidationError("You can only perform these actions to a top-level post")
        if action in (TOGGLE_ACCEPT, MOVE_TO_COMMENT) and self.post.type != Post.ANSWER:
            raise forms.ValidationError("You can only perform these actions to an answer.")
        if action == MOVE_TO_ANSWER and self.post.type != Post.COMMENT:
            raise forms.ValidationError("You can only perform these actions to a comment.")

        if action == DUPLICATE and not dupe:
            raise forms.ValidationError("Unable to close duplicate. Please fill in the post numbers")

        if dupe:
            dupe = dupe.replace(",", " ")
            dupes = dupe.split()[:5]
            cleaned_data['dupe'] = dupes

        return cleaned_data




