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

class EntryList(LoginRequiredMixin, StaffuserRequiredMixin, ListView):
    model = EmailEntry
    template_name = "newsletter/entry_list.html"
    context_object_name = 'entries'

    def get_context_data(self, **kwargs):
        context = super(EntryList, self).get_context_data(**kwargs)
        return context


class EditForm(forms.Form):

    title = forms.CharField(label="Post Title", max_length=200, min_length=10,
                            help_text="Descriptive titles promote better answers.")

    content = forms.CharField(widget=forms.Textarea,
                              min_length=80, max_length=15000,
                              label="Enter your post below")

    text = forms.CharField(widget=forms.Textarea,
                              min_length=30, max_length=15000,
                              label="Text mode")

    status = forms.ChoiceField(choices=((EmailEntry.DRAFT,"Preview"), (EmailEntry.PUBLISHED, "Publish")))

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

DEFAULT_TAGS="biostar exporer"


def now():
    return datetime.datetime.utcnow().replace(tzinfo=utc)

def generate_title():
    count = EmailEntry.objects.all().count() + 1
    date = now().strftime("%b %d, %Y")
    return "Biostar Explorer #%s on %s" % (count, date)

def generate_content():
    html = """
    It's starting to get a little warmer! Tons of new releases
    this week, new version of Python and beta release of Django 1.7! Huzzah!

    Share an article with us and if it lands in newsletter get highlighted as a
    contributor in the newsletter!
    """
    text = "simple text content" * 5
    return html, text

class EditEntry(LoginRequiredMixin, StaffuserRequiredMixin, FormView):
    form_class = EditForm
    model = EmailEntry
    template_name = "newsletter/update_entry.html"

    def get(self, request, *args, **kwargs):
        form = self.form_class()
        pk = int(self.kwargs['pk'])

        if pk > 0:
            # disabled for now
            return Http404
        else:
            title = generate_title()
            html, text = generate_content()

        form = self.form_class(initial=dict(content=html, text=text, title=title))
        return render(request, self.template_name, {'form': form})

    def post(self, request,  *args, **kwargs):
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

