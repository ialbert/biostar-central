# Create your views here.
from django.shortcuts import render_to_response
from django.views.generic import TemplateView, DetailView, ListView, FormView, UpdateView
from .models import User
from . import auth
from django import forms
from django.core.urlresolvers import reverse
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Fieldset, Submit, ButtonHolder
from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect
from django.contrib import messages
from django.core.validators import validate_email
from biostar import const
from braces.views import LoginRequiredMixin
from django.contrib.auth import authenticate, login, logout
from django.conf import settings
from biostar.apps import util
import logging, hmac

logger = logging.getLogger(__name__)

class UserEditForm(forms.Form):
    name = forms.CharField()

    email = forms.EmailField()

    location = forms.CharField(required=False,
                               help_text="Location/Institution (optional)")

    website = forms.URLField(required=False, max_length=200,
                             help_text="The URL to your website (optional)")

    scholar = forms.CharField(required=False, max_length=15,
                              help_text="Your Google Scholar ID (optional)")

    my_tags = forms.CharField(max_length=200, required=False,
                              help_text="Use <code>+</code> to add tags. Add a <code>!</code> to remove a tag. Example: <code>galaxy + bed + solid!</code> (optional)")

    message_prefs = forms.ChoiceField(required=True, choices=const.MESSAGING_TYPE_CHOICES, label="Notifications",
                                      help_text="The default setting for notifications when you contribute to a thread")

    info = forms.CharField(widget=forms.Textarea, required=False,
                           help_text="A brief description about yourself (optional)")

    def __init__(self, *args, **kwargs):
        super(UserEditForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.error_text_inline = False
        self.helper.help_text_inline = True

        self.helper.layout = Layout(
            Fieldset(
                'Update your profile',
                'name',
                'email',
                'location',
                'website',
                'scholar',
                'message_prefs',
                Field('my_tags'),
                'info',
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )


class EditUser(LoginRequiredMixin, FormView):
    """
    Edits a user.
    """

    # The template_name attribute must be specified in the calling apps.
    template_name = ""
    form_class = UserEditForm
    user_fields = "name email".split()
    prof_fields = "location website info scholar my_tags message_prefs".split()

    def get(self, request, *args, **kwargs):
        target = User.objects.get(pk=self.kwargs['pk'])

        initial = {}

        for field in self.user_fields:
            initial[field] = getattr(target, field)

        for field in self.prof_fields:
            initial[field] = getattr(target.profile, field)

        form = self.form_class(initial=initial)
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):
        target = User.objects.get(pk=self.kwargs['pk'])
        target = auth.user_permissions(request=request, target=target)

        # The essential authentication step.
        if not target.has_ownership:
            messages.error(request, "Only owners may edit their profiles")
            return HttpResponseRedirect(reverse("home"))

        form = self.form_class(request.POST)
        if form.is_valid():
            f = form.cleaned_data

            if User.objects.filter(email=f['email']).exclude(pk=request.user.id):
                # Changing email to one that already belongs to someone else.
                messages.error(request, "The email that you've entered is already registered to another user!")
                return render(request, self.template_name, {'form': form})

            # Valid data. Save model attributes and redirect.
            for field in self.user_fields:
                setattr(target, field, f[field])

            for field in self.prof_fields:
                setattr(target.profile, field, f[field])

            target.save()
            target.profile.save()
            messages.success(request, "Profile updated")
            return HttpResponseRedirect(self.get_success_url())

        # There is an error in the form.
        return render(request, self.template_name, {'form': form})

    def get_success_url(self):
        return reverse("user-details", kwargs=dict(pk=self.kwargs['pk']))


def test_login(request):
    # Used during debugging external authentication
    response = redirect("/")
    for name, key in settings.EXTERNAL_AUTH:
        email = "foo@bar.com"
        digest = hmac.new(key, email).hexdigest()
        value = "%s:%s" % (email, digest)
        response.set_cookie(name, value)
        messages.info(request, "set cookie %s, %s, %s" % (name, key, value))
    return response


def external_logout(request):
    "This is required to invalidate the external logout cookies"
    logout(request)
    url = settings.EXTERNAL_LOGOUT_URL or 'account_logout'
    response = redirect(url)
    for name, key in settings.EXTERNAL_AUTH:
        response.delete_cookie(name)
    return response


def external_login(request):
    "This is required to allow external login to proceed"
    url = settings.EXTERNAL_LOGIN_URL or 'account_login'
    response = redirect(url)
    return response

# Adding a captcha enabled form
from allauth.account.views import SignupForm, SignupView
from biostar.apps.util.captcha.fields import MathCaptchaField


class CaptchaForm(SignupForm):
    captcha = MathCaptchaField()


class CaptchaView(SignupView):
    form_class = CaptchaForm

    def get_form_class(self):
        # This is to allow tests to override the form class during testing.
        if settings.CAPTCHA:
            return CaptchaForm
        else:
            return SignupForm

