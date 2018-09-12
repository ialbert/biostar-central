
import logging
import mistune
from django import forms

from django.contrib import messages
from snowpenguin.django.recaptcha2.fields import ReCaptchaField
from snowpenguin.django.recaptcha2.widgets import ReCaptchaWidget
from django.contrib.auth.models import User
from django.conf import settings
from pagedown.widgets import PagedownWidget
from .models import Profile
from . import auth


logger = logging.getLogger("engine")


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

        # TODO: needs to change
        user.username = name.split()[0] + str(user.profile.uid)
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


class EditProfile(forms.Form):

    def __init__(self, user,  *args, **kwargs):

        self.user = user

        super(EditProfile, self).__init__(*args, **kwargs)

        self.fields["email"] = forms.CharField(label='Email', max_length=100, initial=self.user.email)
        self.fields["name"] = forms.CharField(label='Name', max_length=100, initial=self.user.profile.name)
        self.fields["username"] = forms.CharField(label="Handler", max_length=100, initial=self.user.username)
        self.fields["location"] = forms.CharField(label="Location", max_length=100, initial=self.user.profile.location,
                                                  required=False)
        self.fields["website"] = forms.URLField(label="Website", max_length=225, initial=self.user.profile.website,
                                                required=False)
        self.fields["twitter"] = forms.CharField(label="Twitter Id", max_length=100, initial=self.user.profile.twitter,
                                                 required=False)
        self.fields["scholar"] = forms.CharField(label="Scholar", max_length=100, initial=self.user.profile.scholar,
                                                 required=False)
        self.fields["text"] = forms.CharField(widget=PagedownWidget(template="widgets/pagedown.html"),
                               min_length=2, max_length=5000, initial=self.user.profile.text, required=False)
        self.fields["my_tags"] = forms.CharField(max_length=100, required=False, initial=self.user.profile.my_tags,
                                  help_text="""Post with tags listed here will show up in the My Tags tab. 
                                  Use a comma to separate tags. 
                                  Add a <code>!</code> to remove a tag. Example: <code>galaxy, bed, solid!</code> (optional)
                                  """
                                  )
        self.fields["digest_prefs"] = forms.ChoiceField(required=True, choices=Profile.DIGEST_CHOICES, label="Email Digest",
                                         help_text="""(This feature is not working yet!). 
                                         Sets the frequence of digest emails. 
                                         A digest email is a summary of events on the site.""",
                                         initial=self.user.profile.digest_prefs
                                         )

        self.fields["message_prefs"] = forms.ChoiceField(required=True, choices=Profile.MESSAGING_TYPE_CHOICES, label="Notifications",
                                          help_text="""Default mode  sends you an email 
                                          if you receive anwers to questions that you've posted.""",
                                          initial=self.user.profile.message_prefs)

    def save(self, request):

        email = self.cleaned_data['email']
        email_verified = self.user.profile.email_verified

        if self.user.email != email:
            messages.info(request, "Email has to be reverified since it has changed")
            email_verified = False

        self.user.email = self.cleaned_data['email']


        self.user.username = self.cleaned_data["username"]
        self.user.save()

        Profile.objects.filter(user=self.user).update(name=self.cleaned_data['name'],
                                                      location=self.cleaned_data['location'],
                                                      website=self.cleaned_data['website'],
                                                      twitter=self.cleaned_data['twitter'],
                                                      scholar=self.cleaned_data['scholar'],
                                                      text=self.cleaned_data["text"],
                                                      my_tags=self.cleaned_data["my_tags"],
                                                      digest_prefs=self.cleaned_data["digest_prefs"],
                                                      message_prefs=self.cleaned_data["message_prefs"],
                                                      html=mistune.markdown(self.cleaned_data["text"]),
                                                      email_verified=email_verified)
        return self.user

    def clean_email(self):

        data = self.cleaned_data['email']
        email = User.objects.exclude(pk=self.user.pk).filter(email=data)

        if email.exists():
            raise forms.ValidationError("This email is already being used.")

        return data

    def clean_username(self):

        data = self.cleaned_data['username']
        username = User.objects.exclude(pk=self.user.pk).filter(username=data)

        if len(data.split()) > 1:
            raise forms.ValidationError("No spaces allowed in username/handlers.")
        if username.exists():
            raise forms.ValidationError("This handler is already being used.")

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





