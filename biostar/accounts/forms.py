from django import forms
from snowpenguin.django.recaptcha2.fields import ReCaptchaField
from snowpenguin.django.recaptcha2.widgets import ReCaptchaWidget
from django.contrib.auth.models import User
from django.conf import settings
from pagedown.widgets import PagedownWidget
from .models import Profile


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


class SignUpWithCaptcha(SignUpForm):

    def __init__(self, *args, **kwargs):
        super(SignUpWithCaptcha, self).__init__(*args, **kwargs)

        if settings.RECAPTCHA_PRIVATE_KEY:
            self.fields["captcha"] = ReCaptchaField(widget=ReCaptchaWidget())


class LogoutForm(forms.Form):
    pass


class EditProfile(forms.Form):

    email = forms.CharField(label='Email', max_length=100)
    name = forms.CharField(label='Name', max_length=100)
    username = forms.CharField(label="Handler", max_length=100)
    location = forms.CharField(label="Location", max_length=100, required=False)
    website = forms.URLField(label="Website", max_length=225, required=False)
    twitter = forms.CharField(label="Twitter Id", max_length=100, required=False)
    scholar = forms.CharField(label="Scholar", max_length=100, required=False)
    text = forms.CharField(widget=PagedownWidget(template="widgets/pagedown.html"),
                           min_length=2, max_length=5000)
    my_tags = forms.CharField(max_length=100, required=False,
                              help_text="""Post with tags listed here will show up in the My Tags tab. 
                              Use a comma to separate tags. 
                              Add a <code>!</code> to remove a tag. Example: <code>galaxy, bed, solid!</code> (optional)
                              """
                              )
    digest_prefs = forms.ChoiceField(required=True, choices=Profile.DIGEST_CHOICES, label="Email Digest",
                                     help_text="""(This feature is not working yet!). 
                                     Sets the frequence of digest emails. 
                                     A digest email is a summary of events on the site."""
                                     )

    message_prefs = forms.ChoiceField(required=True, choices=Profile.MESSAGING_TYPE_CHOICES, label="Notifications",
                                      help_text="""Default mode  sends you an email 
                                      if you receive anwers to questions that you've posted.""")

    def __init__(self, user,  *args, **kwargs):

        self.user = user

        super(EditProfile, self).__init__(*args, **kwargs)

    def save(self):

        self.user.email = self.cleaned_data['email']
        self.user.profile.name = self.cleaned_data['name']
        self.user.profile.location = self.cleaned_data['location']
        self.user.profile.website = self.cleaned_data['website']

        self.user.save()

        return self.user

    def clean_email(self):

        data = self.cleaned_data['email']
        email = User.objects.exclude(pk=self.user.pk).filter(email=data)

        if email.exists():
            raise forms.ValidationError("This email is already being used.")

        return data

    def clean_username(self):

        data = self.cleaned_data['username']
        username = User.objects.exclude(pk=self.user.pk).filter(username=self.data)

        if username.exists():
            raise forms.ValidationError("This username is already being used.")

        return data


class LoginForm(forms.Form):

    email = forms.CharField(label='Email', max_length=100)
    password = forms.CharField(label='Password', max_length=100,
                               widget=forms.PasswordInput)


class UserModerate(forms.Form):

    CHOICES = [
        (Profile.NEW, "Reinstate as new user"),
        (Profile.TRUSTED, "Reinstate as trusted user"),
        (Profile.BANNED, "Ban user"),
        (Profile.SUSPENDED, "Suspend user")
    ]

    action = forms.ChoiceField(choices=CHOICES, widget=forms.RadioSelect(), label="Select Action")

    def __init__(self, source, target, request, *args, **kwargs):
        self.source = source
        self.target = target
        self.request = request
        super(UserModerate, self).__init__(*args, **kwargs)

    def save(self):
        cleaned_data = self.cleaned_data
        state = cleaned_data["action"]

        Profile.objects.filter(user=self.target).update(state=state)

        return

    def clean(self):
        cleaned_data = super(UserModerate, self).clean()
        if not self.source.profile.is_moderator:
            forms.ValidationError("You need to be a moderator to perform that action")





