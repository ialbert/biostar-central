from django import forms
from .models import Post
from django.core.exceptions import ValidationError
from django.conf import settings

from . import  models, auth
from pagedown.widgets import PagedownWidget
import langdetect

# Share logger with models
logger = models.logger


def english_only(text):


    try:
        text.encode('ascii')
    except Exception:
        raise ValidationError('Title may only contain plain text (ASCII) characters')


def valid_language(text):
    supported_languages = settings.LANGUAGE_DETECTION
    if supported_languages:
        try:
            lang = langdetect.detect(text)
        except Exception as exc:
            logger.error(f"{exc}")
            lang = None

        if lang not in supported_languages:
            raise ValidationError(
                    f'Language "{lang}" is not one of the supported languages {supported_languages}!')


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

    POST_CHOICES = [(Post.QUESTION, "Question"),
                    (Post.JOB, "Job Ad"),
                    (Post.TUTORIAL, "Tutorial"), (Post.TOOL, "Tool"),
                    (Post.FORUM, "Forum"), (Post.NEWS, "News"),
                    (Post.BLOG, "Blog"), (Post.PAGE, "Page")]

    title = forms.CharField(
        label="Post Title",
        max_length=200, min_length=2, validators=[valid_title, english_only],
        help_text="Descriptive titles promote better answers.")

    post_type = forms.ChoiceField(
        label="Post Type",
        choices=POST_CHOICES, help_text="Select a post type: Question, Forum, Job, Blog")

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
                                tag_val=tag_val, author=author)

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
        sub = models.Subscription.get_sub(post=self.post, user=self.user).first()

        if sub:
            sub.type = sub_type
            sub.save()
        else:
            sub = auth.create_sub(post=self.post, user=self.user, sub_type=sub_type)

        return sub



class PostShortForm(forms.Form):
    content = forms.CharField(widget=PagedownWidget(template="widgets/pagedown.html"),
                              min_length=2, max_length=5000)

    def save(self, parent, author, post_type=Post.ANSWER):
        data = self.cleaned_data.get
        answer = auth.create_post(title=parent.title,
                                  parent=parent,
                                  author=author,
                                  content=data("content"),
                                  post_type=post_type
                                  )
        return answer







