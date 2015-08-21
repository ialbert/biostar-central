# Create your views here.
from django.shortcuts import render_to_response
from django.views.generic import TemplateView, DetailView, ListView, FormView, UpdateView
from .models import User, Profile
from . import auth
from django import forms
from django.core.urlresolvers import reverse
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Fieldset, Submit, ButtonHolder, Div
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
    name = forms.CharField(help_text="The name displayed on the site (required)")

    email = forms.EmailField(help_text="Your email, it will not be visible to other users (required)")

    location = forms.CharField(required=False,
                               help_text="Country/City/Institution (recommended)")

    website = forms.URLField(required=False, max_length=100,
                             help_text="The URL to your website (optional)")

    twitter_id = forms.CharField(required=False, max_length=15,
                                 help_text="Your twitter id (optional)")

    scholar = forms.CharField(required=False, max_length=15,
                              help_text="Your Google Scholar ID (optional)")

    my_tags = forms.CharField(max_length=100, required=False,
                              help_text="Post with tags listed here will show up in the My Tags tab. Use a comma to separate tags. Add a <code>!</code> to remove a tag. Example: <code>galaxy, bed, solid!</code> (optional)")

    watched_tags = forms.CharField(max_length=100, required=False,
                                   help_text="Get email when a post matching the tag is posted. Example: <code>minia, bedops, breakdancer, music</code>.")

    digest_prefs = forms.ChoiceField(required=True, choices=Profile.DIGEST_CHOICES, label="Email Digest",
                                     help_text="(This feature is not working yet!). Sets the frequence of digest emails. A digest email is a summary of events on the site.")

    message_prefs = forms.ChoiceField(required=True, choices=const.MESSAGING_TYPE_CHOICES, label="Notifications",
                                      help_text="Default mode  sends you an email if you receive anwers to questions that you've posted.")

    info = forms.CharField(widget=forms.Textarea, required=False, label="Add some information about yourself",
                           help_text="A brief description about yourself (recommended)")

    def __init__(self, *args, **kwargs):
        super(UserEditForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.error_text_inline = False
        self.helper.help_text_inline = True
        self.helper.layout = Layout(
            Fieldset(
                'Update your profile',
                Div(
                    Div('name', ),
                    Div('email', ),
                    Div('location'),
                    Div('website'),
                    css_class="col-md-offset-3 col-md-6",
                ),
                Div(
                    Div('twitter_id', css_class="col-md-6"),
                    Div('scholar', css_class="col-md-6"),
                    Div('digest_prefs', css_class="col-md-6"),
                    Div('message_prefs', css_class="col-md-6"),
                    css_class="col-md-12",
                ),
                Div(
                    Div('my_tags'),
                    Div('watched_tags'),
                    Div('info'),
                    ButtonHolder(
                        Submit('submit', 'Submit')
                    ),
                    css_class="col-md-12",
                ),
            ),

        )


class EditUser(LoginRequiredMixin, FormView):
    """
    Edits a user.
    """

    # The template_name attribute must be specified in the calling apps.
    template_name = ""
    form_class = UserEditForm
    user_fields = "name email".split()
    prof_fields = "location website info scholar my_tags watched_tags twitter_id message_prefs digest_prefs".split()

    def get(self, request, *args, **kwargs):
        target = User.objects.get(pk=self.kwargs['pk'])
        target = auth.user_permissions(request=request, target=target)
        if not target.has_ownership:
            messages.error(request, "Only owners may edit their profiles")
            return HttpResponseRedirect(reverse("home"))

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
        profile = target.profile

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
                setattr(profile, field, f[field])

            target.save()
            profile.add_tags(profile.watched_tags)
            profile.save()
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


class DigestForm(forms.Form):
    digest_prefs = forms.ChoiceField(required=True, choices=Profile.DIGEST_CHOICES, label="Emails",
                                     help_text="Sets the frequence of digest emails")


    def __init__(self, *args, **kwargs):
        super(DigestForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.error_text_inline = False
        self.helper.help_text_inline = True
        self.helper.layout = Layout(
            Fieldset(
                'Update your digest settings',
                Div(
                    Div('digest_prefs', ),
                    ButtonHolder(
                        Submit('submit', 'Submit')
                    ),
                    css_class="col-md-6",
                ),
            ),
        )

def unsubscribe(request, uuid):

    user = User.objects.filter(profile__uuid=uuid)
    if not user:
        messages.error(request, 'Invalid user identifier.')
    else:
        user[0].profile.digest_prefs = Profile.NO_DIGEST
        messages.success(request, 'User unsubscribed.')

    response = redirect(reverse("home"))
    return response

class DigestManager(LoginRequiredMixin, FormView):
    "Handle user digest subscriptions"

    template_name = "digest/manager.html"
    form_class = DigestForm

    def get(self, request, *args, **kwargs):
        form = self.form_class()
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):
        form = self.form_class(request.POST)
        if form.is_valid():
            f = form.cleaned_data
            user = request.user
            user.profile.digest_prefs = int(form.cleaned_data['digest_prefs'])
            user.profile.save()
            messages.success(request, 'Email frequency has been set to %s' % user.profile.get_digest_prefs_display())

        return render(request, self.template_name, {'form': form})


def external_login(request):
    "This is required to allow external login to proceed"
    url = settings.EXTERNAL_LOGIN_URL or 'account_login'
    response = redirect(url)
    return response

# Adding a captcha enabled form
from allauth.account.views import SignupForm, SignupView
from biostar.apps.util.captcha.fields import MathCaptchaField
from captcha.fields import ReCaptchaField


class CaptchaForm(SignupForm):
    captcha = ReCaptchaField()


class CaptchaView(SignupView):
    form_class = CaptchaForm

    def get_form_class(self):

        # This is to allow tests to override the form class during testing.
        if settings.RECAPTCHA_PRIVATE_KEY:
            return CaptchaForm
        else:
            return SignupForm
