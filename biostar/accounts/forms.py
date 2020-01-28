
import logging
from django import forms

from django.contrib import messages
from snowpenguin.django.recaptcha2.fields import ReCaptchaField
from snowpenguin.django.recaptcha2.widgets import ReCaptchaWidget
from django.contrib.auth.models import User
from django.conf import settings
from .models import Profile
from . import auth, util


logger = logging.getLogger("engine")

MAX_TAGS = 50


class SignUpForm(forms.Form):

    password1 = forms.CharField(
        label="Password",
        strip=False,
        widget=forms.PasswordInput,
        max_length=254,
        min_length=2,
    )

    password2 = forms.CharField(
        label="Password confirmation",
        widget=forms.PasswordInput,
        strip=False,
        help_text="Enter the same password as before, for verification.",
    )

    email = forms.CharField(
        label="Email",
        strip=False,
        widget=forms.TextInput,
        max_length=254,
        min_length=2,
    )

    def clean_password2(self):

        password1 = self.cleaned_data.get("password1")
        password2 = self.cleaned_data.get("password2")
        if password1 and password2 and password1 != password2:
            raise forms.ValidationError("Passwords given do not match.")
        return password2

    def clean_email(self):

        data = self.cleaned_data['email']
        if User.objects.filter(email=data).exists():
            raise forms.ValidationError("This email is already being used.")
        return data

    def save(self):

        email = self.cleaned_data.get('email')
        password = self.cleaned_data.get('password1')
        name = email.split("@")[0]
        user = User.objects.create(email=email, first_name=name)
        user.set_password(password)
        user.save()

        # Send
        auth.send_verification_email(user=user)
        logger.info(f"Signed up user.id={user.id}, user.email={user.email}")

        return user


class SignUpWithCaptcha(SignUpForm):

    def __init__(self, *args, **kwargs):
        super(SignUpWithCaptcha, self).__init__(*args, **kwargs)

        if settings.RECAPTCHA_PRIVATE_KEY:
            self.fields["captcha"] = ReCaptchaField(widget=ReCaptchaWidget())


class LogoutForm(forms.Form):
    pass


def validate_tags(tags):
    my_tags = tags.split(',')
    if len(my_tags) > MAX_TAGS:
        return forms.ValidationError("Maximum number of tags reached.")
    return tags


class EditProfile(forms.Form):
    name = forms.CharField(label='Name', max_length=100, required=True)
    email = forms.CharField(label='Email', max_length=100, required=True)
    username = forms.CharField(label="Handler", max_length=100, required=True)
    location = forms.CharField(label="Location", max_length=100, required=False)
    website = forms.URLField(label="Website", max_length=225, required=False)
    twitter = forms.CharField(label="Twitter Id", max_length=100, required=False)
    scholar = forms.CharField(label="Scholar", max_length=100, required=False)

    text = forms.CharField(widget=forms.Textarea(),min_length=2, max_length=5000, required=False,
                           help_text="Extra information about you to personalize your profile.")

    message_prefs = forms.ChoiceField(required=True, label="Notifications", choices=Profile.MESSAGING_TYPE_CHOICES,
                                      widget=forms.Select(attrs={'class': "ui dropdown"}),
                                      help_text="""Default mode sends notifications using local messages.""")
    my_tags = forms.CharField(label="My tags", max_length=500, required=False,
                              help_text="""
                              Add a tag by typing a word then adding a comma or press ENTER or SPACE.
                              """, widget=forms.HiddenInput())
    watched_tags = forms.CharField(label="Watched tags", max_length=50, required=False,
                                   help_text="""
                              Add a tag by typing a word then adding a comma or press ENTER or SPACE.
                              """, widget=forms.HiddenInput())

    def __init__(self, user,  *args, **kwargs):

        self.user = user

        super(EditProfile, self).__init__(*args, **kwargs)

    def clean_username(self):

        data = self.cleaned_data['username']
        username = User.objects.exclude(pk=self.user.pk).filter(username=data)

        if len(data.split()) > 1:
            raise forms.ValidationError("No spaces allowed in username/handlers.")
        if username.exists():
            raise forms.ValidationError("This handler is already being used.")

        return data

    def clean_email(self):
        cleaned_data = self.cleaned_data['email']
        email = User.objects.filter(email=cleaned_data).exclude(pk=self.user.pk).first()

        if email:
            raise forms.ValidationError("Email already exists.")

        if self.user.is_superuser and cleaned_data != self.user.email:
            raise forms.ValidationError("Admins are required to change emails using the Django Admin Interface.")

        return cleaned_data

    def clean_my_tags(self):
        my_tags = self.cleaned_data['my_tags']
        my_tags = ','.join(list(set(my_tags.split(","))))
        return validate_tags(tags=my_tags)

    def clean_watched_tags(self):
        watched_tags = self.cleaned_data['watched_tags']
        watched_tags = ','.join(list(set(watched_tags.split(","))))
        return validate_tags(tags=watched_tags)


class LoginForm(forms.Form):

    email = forms.CharField(label='Email', max_length=100)
    password = forms.CharField(label='Password', max_length=100,
                               widget=forms.PasswordInput)


class UserModerate(forms.Form):

    CHOICES = [
        (Profile.SPAMMER, 'Report as spammer'),
        (Profile.BANNED, "Ban user"),
        (Profile.SUSPENDED, "Suspend user"),
        (Profile.NEW, "Reinstate as new user"),
        (Profile.TRUSTED, "Reinstate as trusted user")
    ]

    action = forms.IntegerField(widget=forms.RadioSelect(choices=CHOICES), required=False, label="Select Action")

    def __init__(self, source, target, request, *args, **kwargs):
        self.source = source
        self.target = target
        self.request = request
        super(UserModerate, self).__init__(*args, **kwargs)

    def clean(self):
        cleaned_data = super(UserModerate, self).clean()
        action = cleaned_data['action']

        if not self.source.profile.is_moderator:
            raise forms.ValidationError("You need to be a moderator to perform that action")

        if action == Profile.BANNED and not self.source.is_superuser:
            raise forms.ValidationError("You need to be an admin to ban users.")

        if self.target.profile.is_moderator and not self.source.is_superuser:
            raise forms.ValidationError("You need to be an admin to moderator other moderators.")

        if self.target == self.source:
            raise forms.ValidationError("You can not moderate yourself.")

