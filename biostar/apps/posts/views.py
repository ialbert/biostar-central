# Create your views here.
from django.shortcuts import render_to_response
from django.views.generic import TemplateView, DetailView, ListView, FormView, UpdateView
from .models import Post, Torrent
from django import forms
from django.core.urlresolvers import reverse
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Fieldset, Div, Submit, ButtonHolder
from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.contrib import messages
from . import auth
from braces.views import LoginRequiredMixin
from datetime import datetime
from django.utils.timezone import utc
from django.core.exceptions import ObjectDoesNotExist
from django.conf import settings
from biostar.const import OrderedDict
from django.core.exceptions import ValidationError
from django.contrib.auth.decorators import login_required
from django.views.decorators.csrf import csrf_exempt
import logging
from biostar.apps.tracker.utils import torrent_get_hash

logger = logging.getLogger(__name__)


def valid_title(text):
    "Validates form input for tags"
    text = text.strip()
    if not text:
        raise ValidationError('Please enter a title')

    if len(text) < 10:
        raise ValidationError('The title is too short')

    words = text.split(" ")
    if len(words) < 3:
        raise ValidationError('More than two words please.')


def valid_tag(text):
    "Validates form input for tags"
    text = text.strip()
    if not text:
        raise ValidationError('Please enter at least one tag')
    if len(text) > 50:
        raise ValidationError('Tag line is too long (50 characters max)')
    words = text.split(",")
    if len(words) > 5:
        raise ValidationError('You have too many tags (5 allowed)')


class DataForm(forms.Form):
    data = forms.FileField(
        label="Optional: Share Data",
        required=False,
        help_text='File must be a torrent file! Read the <a href="/info/data" target="blank">Data Sharing</a> for how this works'
    )

    def clean_data(self):
        MAX_MB = 1
        MAX_SIZE = MAX_MB * 1024 * 1024
        data = self.cleaned_data['data']
        if not data:
            return data

        if data._size > MAX_SIZE:
            raise forms.ValidationError('Filesize must be under %s Mb' % MAX_MB)

        try:
            value = data.read()
            info_hash = torrent_get_hash(value)
        except Exception, exc:
            logger.error(exc)
            raise forms.ValidationError('The file does not appear to be a torrent file')

        return data


class LongForm(DataForm):
    FIELDS = "title content post_type tag_val".split()

    POST_CHOICES = [(Post.QUESTION, "Question"),
                    (Post.JOB, "Job Ad"),
                    (Post.DATA, "Data"),
                    (Post.TUTORIAL, "Tutorial"), (Post.TOOL, "Tool"),
                    (Post.FORUM, "Forum"), (Post.NEWS, "News"),
                    (Post.BLOG, "Blog"), ]

    title = forms.CharField(
        label="Post Title",
        max_length=200, min_length=10, validators=[valid_title],
        help_text="Descriptive titles promote better answers.")

    post_type = forms.ChoiceField(
        label="Post Type",
        choices=POST_CHOICES, help_text="Select a post type: Question, Forum, Job, Blog")

    tag_val = forms.CharField(
        label="Post Tags",
        required=True, validators=[valid_tag],
        help_text="Choose one or more tags to match the topic. To create a new tag just type it in and press ENTER.",
    )

    content = forms.CharField(widget=forms.Textarea,
                              min_length=80, max_length=15000,
                              label="Enter your post below")

    remove = forms.BooleanField(
        label="Remove the data from this post", initial=False, required=False,
        help_text="Checking this box will remove the data. Uploading a new file will replace the data",
    )

    POST_KEY = "post"

    def __init__(self, *args, **kwargs):
        # This is sent only when editing a post
        post = kwargs.get(self.POST_KEY)
        if post:
            kwargs.pop(self.POST_KEY)

        super(LongForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_class = "post-form"
        self.helper.layout = Layout(
            Fieldset(
                'Post Form',
                Field('title'),
                Field('post_type'),
                Field('tag_val'),
                Field('content'),
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )

        fieldset = self.helper.layout[0]

        #if post and post.torrent:
        #    fieldset.append(
        #        Field('remove'),
        #    )

        # add data upload
        fieldset.append(
            Field('data'),
        )


class ShortForm(DataForm):
    FIELDS = ["content"]

    content = forms.CharField(widget=forms.Textarea, min_length=20)

    data = forms.FileField(
        label="Optional: Share Data",
        required=False,
        help_text='File must be a torrent file! Read the <a href="/info/data" target="blank">Data Sharing</a> for how this works'
    )

    POST_KEY = "post"

    def __init__(self, *args, **kwargs):
        post = kwargs.get(self.POST_KEY)
        if post:
            kwargs.pop(self.POST_KEY)
        super(ShortForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.layout = Layout(
            Fieldset(
                'Post',
                Field('content'),
                Field('data'),
            ),

            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )


def parse_tags(category, tag_val):
    pass


@login_required
@csrf_exempt
def external_post_handler(request):
    "This is used to pre-populate a new form submission"
    import hmac

    user = request.user
    home = reverse("home")
    name = request.REQUEST.get("name")

    if not name:
        messages.error(request, "Incorrect request. The name parameter is missing")
        return HttpResponseRedirect(home)

    try:
        secret = dict(settings.EXTERNAL_AUTH).get(name)
    except Exception, exc:
        logger.error(exc)
        messages.error(request, "Incorrect EXTERNAL_AUTH settings, internal exception")
        return HttpResponseRedirect(home)

    if not secret:
        messages.error(request, "Incorrect EXTERNAL_AUTH, no KEY found for this name")
        return HttpResponseRedirect(home)

    content = request.REQUEST.get("content")
    submit = request.REQUEST.get("action")
    digest1 = request.REQUEST.get("digest")
    digest2 = hmac.new(secret, content).hexdigest()

    if digest1 != digest2:
        messages.error(request, "digests does not match")
        return HttpResponseRedirect(home)

    # auto submit the post
    if submit:
        post = Post(author=user, type=Post.QUESTION)
        for field in settings.EXTERNAL_SESSION_FIELDS:
            setattr(post, field, request.REQUEST.get(field, ''))
        post.save()
        post.add_tags(post.tag_val)
        return HttpResponseRedirect(reverse("post-details", kwargs=dict(pk=post.id)))

    # pre populate the form
    sess = request.session
    sess[settings.EXTERNAL_SESSION_KEY] = dict()
    for field in settings.EXTERNAL_SESSION_FIELDS:
        sess[settings.EXTERNAL_SESSION_KEY][field] = request.REQUEST.get(field, '')

    return HttpResponseRedirect(reverse("new-post"))


def add_data(post, data):
    if not data:
        return

    # rewind the incoming file
    data.seek(0)
    content = data.read()
    name = data.name

    info_hash = torrent_get_hash(content)
    torrent = Torrent.objects.create(
        post=post,
        info_hash=info_hash,
        name=name, content=content,
    )
    post.data_count = post.torrent_set.all().count()
    post.save()


class NewPost(LoginRequiredMixin, FormView):
    form_class = LongForm
    template_name = "post_edit.html"

    def get(self, request, *args, **kwargs):
        initial = dict()

        # Attempt to prefill from GET parameters
        for key in "title tag_val content".split():
            value = request.GET.get(key)
            if value:
                initial[key] = value

        # Attempt to prefill from external session
        sess = request.session
        if settings.EXTERNAL_SESSION_KEY in sess:
            for field in settings.EXTERNAL_SESSION_FIELDS:
                initial[field] = sess[settings.EXTERNAL_SESSION_KEY].get(field)
            del sess[settings.EXTERNAL_SESSION_KEY]

        form = self.form_class(initial=initial)
        return render(request, self.template_name, {'form': form})


    def post(self, request, *args, **kwargs):

        # Validating the form.
        form = self.form_class(request.POST, request.FILES)

        if not form.is_valid():
            return render(request, self.template_name, {'form': form})

        # Valid forms start here.
        get = form.cleaned_data.get

        title = get('title')
        content = get('content')
        post_type = int(get('post_type'))
        tag_val = get('tag_val')

        # This the the data set
        data = get('data')

        post = Post(
            title=title, content=content, tag_val=tag_val,
            author=request.user, type=post_type,
        )
        post.save()

        add_data(post=post, data=data)

        # Triggers a new post save.
        post.add_tags(post.tag_val)

        messages.success(request, "%s created" % post.get_type_display())
        return HttpResponseRedirect(post.get_absolute_url())


class NewAnswer(LoginRequiredMixin, FormView):
    """
    Creates a new post.
    """
    form_class = ShortForm
    template_name = "post_edit.html"
    type_map = dict(answer=Post.ANSWER, comment=Post.COMMENT)
    post_type = None

    def get(self, request, *args, **kwargs):
        initial = {}

        # The parent id.
        pid = int(self.kwargs['pid'])
        #form_class = ShortForm if pid else LongForm
        form = self.form_class(initial=initial)
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):

        pid = int(self.kwargs['pid'])

        # Find the parent.
        try:
            parent = Post.objects.get(pk=pid)
        except ObjectDoesNotExist, exc:
            messages.error(request, "The post does not exist. Perhaps it was deleted")
            HttpResponseRedirect("/")

        # Validating the form.
        form = self.form_class(request.POST)
        if not form.is_valid():
            return render(request, self.template_name, {'form': form})

        # Valid forms start here.
        clean = form.cleaned_data.get

        # Figure out the right type for this new post
        post_type = self.type_map.get(self.post_type)
        # Create a new post.
        post = Post(
            title=parent.title, content=clean('content'), author=request.user, type=post_type,
            parent=parent,
        )

        messages.success(request, "%s created" % post.get_type_display())
        post.save()

        return HttpResponseRedirect(post.get_absolute_url())


class EditPost(LoginRequiredMixin, FormView):
    """
    Edits an existing post.
    """

    # The template_name attribute must be specified in the calling apps.
    template_name = "post_edit.html"
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

        initial = dict(title=post.title, content=post.content, post_type=post.type, tag_val=post.tag_val)

        form_class = LongForm if post.is_toplevel else ShortForm
        form = form_class(initial=initial, post=post)
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):

        pk = int(self.kwargs['pk'])
        post = Post.objects.get(pk=pk)
        post = auth.post_permissions(request=request, post=post)

        # For historical reasons we had posts with iframes
        # these cannot be edited because the content would be lost in the front end
        if "<iframe" in post.content:
            messages.error(request, "This post is not editable because of an iframe! Contact if you must edit it")
            return HttpResponseRedirect(post.get_absolute_url())

        # Check and exit if not a valid edit.
        if not post.is_editable:
            messages.error(request, "This user may not modify the post")
            return HttpResponseRedirect(post.get_absolute_url())

        # Posts with a parent are not toplevel
        form_class = LongForm if post.is_toplevel else ShortForm

        form = form_class(request.POST, request.FILES)

        if not form.is_valid():
            # Invalid form submission.
            return render(request, self.template_name, {'form': form})

        # Valid forms start here.
        clean = form.cleaned_data

        # Set the form attributes.
        for field in form_class.FIELDS:
            setattr(post, field, clean[field])

        # TODO: fix this oversight!
        post.type = int(clean.get('post_type', post.type))

        # This is needed to validate some fields.
        post.save()

        # Obtain dataset
        data = clean.get('data')

        # add the data if
        add_data(post=post, data=data)

        print post.torrent_set.all()

        if post.is_toplevel:
            post.add_tags(post.tag_val)

        # Update the last editing user.
        post.lastedit_user = request.user
        post.lastedit_date = datetime.utcnow().replace(tzinfo=utc)
        post.save()
        messages.success(request, "Post updated")

        return HttpResponseRedirect(post.get_absolute_url())

    def get_success_url(self):
        return reverse("user_details", kwargs=dict(pk=self.kwargs['pk']))

