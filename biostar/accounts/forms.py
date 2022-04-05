import logging

import bleach
import mistune
import re
from bleach.callbacks import nofollow
from django import forms
from django.conf import settings
from django.contrib.auth.models import User
from django.core.validators import FileExtensionValidator
from django.template.defaultfilters import slugify
from snowpenguin.django.recaptcha2.fields import ReCaptchaField
from snowpenguin.django.recaptcha2.widgets import ReCaptchaWidget

from .models import Profile, UserImage

logger = logging.getLogger("engine")

MAX_TAGS = 50
IMG_EXTENTIONS = ['jpg',
                  'jpeg',
                  'png',
                  'webp'
                  ]


def check_size(fobj, maxsize=0.3, field=None):
    # maxsize in megabytes!
    error_msg = ''
    try:
        if fobj and fobj.size > maxsize * 1024 * 1024.0:
            curr_size = fobj.size / 1024 / 1024.0
            prefix = f'{field} field : '.capitalize() if field else ''
            error_msg = prefix + f"file too large, {curr_size:0.1f}MB should be < {maxsize:0.1f}MB"

    except Exception as exc:
        error_msg = f"File size validation error: {exc}"

    if error_msg:
        raise forms.ValidationError(error_msg)

    return fobj


def valid_tag(text):
    "Validates form input for tags"

    tag_val = text.replace(',', ' ').split()

    if len(tag_val) > MAX_TAGS:
        return forms.ValidationError("Maximum number of tags reached.")

    if settings.STRICT_TAGS:
        pattern = r'^[A-Za-z0-9-._]+$'
        for tag in tag_val:
            match = re.match(pattern, tag)
            if not match:
                raise forms.ValidationError(f'Invalid characters in tag: {tag}')


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

        logger.info(f"Signed up user.id={user.id}, user.email={user.email}")

        return user


class SignUpWithCaptcha(SignUpForm):

    def __init__(self, *args, **kwargs):
        super(SignUpWithCaptcha, self).__init__(*args, **kwargs)

        if settings.RECAPTCHA_PRIVATE_KEY:
            self.fields["captcha"] = ReCaptchaField(widget=ReCaptchaWidget())


class LogoutForm(forms.Form):
    pass


def markdown(text):
    # Add admin urls.
    html = mistune.markdown(text)

    html = bleach.linkify(text=html, callbacks=[nofollow], skip_tags=['pre', 'code'])

    return html

def get_tags_widget(attrs=None):
    attrs = attrs or {}

    if settings.DROPDOWN_TAGS:
        return forms.HiddenInput()
    else:
        return forms.TextInput(attrs=attrs)


class EditProfile(forms.Form):

    def __init__(self, user, *args, **kwargs):

        self.user = user

        super(EditProfile, self).__init__(*args, **kwargs)

        self.fields['name'] = forms.CharField(label='Name', max_length=100, required=True,
                                              initial=self.user.profile.name)
        self.fields['email'] = forms.CharField(label='Email', max_length=100, required=True,
                                               initial=self.user.email)
        self.fields['handle'] = forms.CharField(label="Handler", max_length=100, required=True,
                                                  initial=self.user.profile.handle)
        self.fields['location'] = forms.CharField(label="Location", max_length=100, required=False,
                                                  initial=self.user.profile.location)
        self.fields['website'] = forms.URLField(label="Website", max_length=225, required=False,
                                                initial=self.user.profile.website)
        self.fields['twitter'] = forms.CharField(label="Twitter Id", max_length=100, required=False,
                                                 initial=self.user.profile.twitter)
        self.fields['scholar'] = forms.CharField(label="Scholar", max_length=100, required=False,
                                                 initial=self.user.profile.scholar)

        self.fields['user_icon'] = forms.ChoiceField(required=False, label="User icon",
                                                     choices=Profile.USER_ICON_CHOICES,
                                                     widget=forms.Select(attrs={'class': "ui dropdown"}),
                                                     initial=self.user.profile.user_icon,
                                                     help_text="User icon type")

        self.fields['text'] = forms.CharField(widget=forms.Textarea(attrs={'rows': 20}),
                                              min_length=2, max_length=5000, required=False,
                                              help_text="Information about you (markdown)",

                                              initial=self.user.profile.text)

        self.fields['message_prefs'] = forms.ChoiceField(required=True, label="Notifications",
                                                         choices=Profile.MESSAGING_TYPE_CHOICES,
                                                         widget=forms.Select(attrs={'class': "ui dropdown"}),
                                                         initial=self.user.profile.message_prefs,
                                                         help_text="Default mode sends notifications using local messages.")

        self.fields['digest_prefs'] = forms.ChoiceField(required=True, label="Digest options",
                                                        choices=Profile.DIGEST_CHOICES,
                                                        widget=forms.Select(attrs={'class': "ui dropdown"}),
                                                        initial=self.user.profile.digest_prefs,
                                                        help_text="Digest are sent through the email provided.")

        self.fields['my_tags'] = forms.CharField(label="My tags", max_length=500, required=False,
                                                 initial=self.user.profile.my_tags, validators=[valid_tag],
                                                 help_text="""
                                  Add a tag by typing a word then adding a comma or press ENTER or SPACE.
                                  """, widget=get_tags_widget())
        self.fields['watched_tags'] = forms.CharField(label="Watched tags", max_length=500, required=False,
                                                      help_text="""
                                  Add a tag by typing a word then adding a comma or press ENTER or SPACE.
                                  """,initial=self.user.profile.watched_tags, validators=[valid_tag],
                                  widget=get_tags_widget())

    def clean_handle(self):

        data = self.cleaned_data['handle']
        data = slugify(data)
        handle = Profile.objects.filter(handle=data).exclude(user=self.user)

        if handle.exists():
            raise forms.ValidationError("This handle is already being used.")

        return data

    def clean_email(self):
        cleaned_data = self.cleaned_data['email']
        email = User.objects.filter(email=cleaned_data).exclude(pk=self.user.pk).first()

        if email:
            raise forms.ValidationError("Email already exists.")

        return cleaned_data

    def clean_my_tags(self):
        my_tags = self.cleaned_data["my_tags"]
        my_tags = self.tag_cleaner(tags=my_tags)
        return my_tags

    def clean_watched_tags(self):
        watched_tags = self.cleaned_data["watched_tags"]
        watched_tags = self.tag_cleaner(tags=watched_tags)
        return watched_tags

    def tag_cleaner(self, tags):

        if settings.STRICT_TAGS:
            tags = tags.replace(',', ' ').split()
        else:
            tags = tags.split(',')

        tags = set(tags)
        tags = ",".join(tags)

        return tags

    def save(self):
        email = self.cleaned_data['email']
        html = markdown(self.cleaned_data["text"])

        # Update usernames and email
        User.objects.filter(pk=self.user.pk).update(email=email)

        # Change email verification status if email changes.
        verified = False if email != self.user.email else self.user.profile.email_verified

        # Update profile attributes
        Profile.objects.filter(user=self.user).update(
            html=html,
            name=self.cleaned_data['name'],
            watched_tags=self.cleaned_data['watched_tags'],
            location=self.cleaned_data['location'],
            handle=self.cleaned_data['handle'],
            website=self.cleaned_data['website'],
            twitter=self.cleaned_data['twitter'],
            scholar=self.cleaned_data['scholar'],
            text=self.cleaned_data["text"],
            email_verified=verified,
            my_tags=self.cleaned_data['my_tags'],
            user_icon=self.cleaned_data['user_icon'],
            message_prefs=self.cleaned_data["message_prefs"],
            digest_prefs=self.cleaned_data['digest_prefs'])
        # Recompute watched tags
        Profile.objects.filter(user=self.user).first().add_watched()


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


class ImageUploadForm(forms.Form):

    def __init__(self, user=None, *args, **kwargs):

        self.user = user
        super(ImageUploadForm, self).__init__(*args, **kwargs)

    image = forms.ImageField(required=True,
                             validators=[FileExtensionValidator(allowed_extensions=IMG_EXTENTIONS)])

    def clean_image(self):

        img = self.cleaned_data['image']

        # Get all images this user has uploaded so far.
        userimg = UserImage.objects.filter(user=self.user)

        # Check for current image size being uploaded.
        check_size(fobj=img, maxsize=settings.MAX_IMAGE_SIZE_MB)

        # Moderators get no limit on images.
        if self.user.is_authenticated and self.user.profile.is_moderator:
            return img

        if userimg.count() >= settings.MAX_IMAGES:
            raise forms.ValidationError("Exceeded the maximum amount of images you can upload.")

        return img

    def save(self):

        # Store the file
        image = self.cleaned_data['image']

        # Create user image object
        userimg = UserImage.objects.create(user=self.user)

        # Save image to database.
        userimg.image.save(image.name, image, save=True)

        return userimg.image.url
