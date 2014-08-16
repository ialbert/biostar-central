__author__ = 'ialbert'

from django.views.generic import TemplateView, DetailView, ListView, FormView, UpdateView, CreateView
from .models import Post, EmailSub, EmailEntry
from braces.views import LoginRequiredMixin, StaffuserRequiredMixin
from django.http import HttpResponseRedirect, Http404
from django import forms
from django.shortcuts import render
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Fieldset, Div, Submit, ButtonHolder
from django.contrib import messages
import datetime
from django.utils.timezone import utc
from biostar.apps.util import html
from django.template import loader, Context, Template, RequestContext
from django.contrib.sites.models import Site
from django.contrib.auth import get_user_model

def text_render(request, name, params):
    "Helper function to render a template"
    tmpl = loader.get_template(name)
    cont = RequestContext(request, params)
    page = tmpl.render(cont)
    return page

class EntryList(LoginRequiredMixin, StaffuserRequiredMixin, ListView):
    model = EmailEntry
    template_name = "explorer/email_list.html"
    context_object_name = 'entries'

    def get_context_data(self, **kwargs):
        context = super(EntryList, self).get_context_data(**kwargs)
        return context


class EditForm(forms.Form):
    title = forms.CharField(label="Email Subject", max_length=200, min_length=10,
                            help_text="Descriptive titles promote better answers.")

    content = forms.CharField(widget=forms.Textarea,
                              min_length=80, max_length=15000,
                              label="HTML Email")

    text = forms.CharField(widget=forms.Textarea,
                           min_length=30, max_length=15000,
                           label="Text Email")

    status = forms.ChoiceField(choices=((EmailEntry.DRAFT, "Preview"), (EmailEntry.PUBLISHED, "Publish")))

    def __init__(self, *args, **kwargs):
        super(EditForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.layout = Layout(
            Fieldset(
                'News Form',
                Field('title'),
                Field('content'),
                Field('text'),
                Field('status'),
            ),
            ButtonHolder(
                Submit('submit', 'Submit'),
            )
        )


DEFAULT_TAGS = "biostar exporer"


def now():
    return datetime.datetime.utcnow().replace(tzinfo=utc)


def generate_title(request):
    count = EmailEntry.objects.all().count() + 1
    date = now().strftime("%b %d, %Y")
    return "Biostar Explorer #%s on %s" % (count, date)


def generate_content(request, days=365*10):
    since1 = now() - datetime.timedelta(days=days)

    most_viewed = Post.objects.filter(creation_date__gt=since1, type__in=Post.TOP_LEVEL).order_by("-view_count")[:5]
    most_active = Post.objects.filter(lastedit_date__gt=since1, type__in=Post.TOP_LEVEL).order_by("-view_count")[:5]


    User = get_user_model()
    since2 = now() - datetime.timedelta(days=365*10)
    new_active_users = User.objects.filter(profile__date_joined__gt=since2, score__gt=0).order_by("-score")[:5]

    # Most active users.
    active_users = User.objects.filter(profile__date_joined__gt=since2, score__gt=0).order_by("-score")[:5]

    protocol = 'https://' if request.is_secure() else 'http://'

    params = dict(
        protocol = protocol,
        site=Site.objects.get_current(),
        most_viewed=most_viewed,
        most_active=most_active,
        new_active_users = new_active_users,
        active_users = active_users,
    )
    rich = text_render(request, "explorer/email_body.html", params)
    text = text_render(request, "explorer/email_body.txt", params)
    return rich, text


class EditEntry(LoginRequiredMixin, StaffuserRequiredMixin, FormView):
    form_class = EditForm
    model = EmailEntry
    template_name = "explorer/email_edit.html"

    def get(self, request, *args, **kwargs):
        form = self.form_class()

        title = generate_title(request)
        rich, text = generate_content(request)

        form = self.form_class(initial=dict(content=rich, text=text, title=title))
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):
        form = self.form_class(request.POST)

        if not form.is_valid():
            messages.error(request, "form is NOT valid")
            return render(request, self.template_name, {'form': form})

        get = form.cleaned_data.get
        title, content, text, status = get('title'), get('content'), get('text'), int(get('status'))

        if status == EmailEntry.PUBLISHED:
            post = Post(title=title, content=content, author=self.request.user)
            post.save()
            post.add_tags(DEFAULT_TAGS)

            # Create the corresponding email entry.
            EmailEntry.objects.create(post=post, status=EmailEntry.PENDING, text=text)

            return HttpResponseRedirect(post.get_absolute_url())

        messages.success(request, "form is valid")
        return render(request, self.template_name, {'form': form})

