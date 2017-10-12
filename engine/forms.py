import hjson as json

from django import forms
from django.utils.translation import gettext_lazy as helpers
from django.contrib.auth.models import User
from .models import Project, Data, Analysis
from . import util
from . import factory

from pagedown.widgets import PagedownWidget


class SignUpForm(forms.ModelForm):
    
    password1 = forms.CharField(
        label=helpers("Password"),
        strip=False,
        widget=forms.PasswordInput,
        max_length=254,
        min_length=2,
    )

    password2 = forms.CharField(
        label=helpers("Password confirmation"),
        widget=forms.PasswordInput,
        strip=False,
        help_text=helpers("Enter the same password as before, for verification."),
    )

    class Meta:

        model = User
        fields = ("email",)

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)


    def clean_password2(self):

        password1 = self.cleaned_data.get("password1")
        password2 = self.cleaned_data.get("password2")
        if password1 and password2 and password1 != password2:
            raise forms.ValidationError(
                helpers("Passwords given do not match."))
        return password2

    def clean_email(self):

        data = self.cleaned_data['email']
        if User.objects.filter(email=data).exists():
            raise forms.ValidationError("This email is already being used.")
        return data

    def cleaned_data(self, *args):
        return self.cleaned_data


class LoginForm(forms.Form):

    email = forms.CharField(label='Email', max_length=100)
    password = forms.CharField(label='Password', max_length=100, 
                              widget=forms.PasswordInput)


class ProjectForm(forms.ModelForm):

    text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))

    class Meta:

        model = Project
        fields = ['title', 'text']


class DataForm(forms.Form):

    title = forms.CharField(max_length=256)
    text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))
    file = forms.FileField(label="Upload data file")
    type = forms.IntegerField(widget=forms.Select(choices=Data.TYPE_CHOICES),
                              initial=Data.FILE)
    def save(self, *args, **kwargs):

        super(DataForm, self).save(*args, **kwargs)


class RunAnalysis(forms.Form):

    def __init__(self, *args, **kwargs):

        analysis = kwargs.pop("analysis")

        super().__init__(*args, **kwargs)

        analysis = json.loads(analysis)
        # Job needs a title
        self.fields["title"] = forms.CharField(max_length=256, initial="Title")

        for field in analysis:

            data = analysis[field]
            display_type = data["display_type"]
            factory.check_display(display_type)

            if data.get("visible") == 1:

                self.fields[field] = factory.TYPE2FUNC[display_type](data)

    def save(self, *args, **kwargs):

        super(RunAnalysis, self).save(*args, **kwargs)


class EditAnalysis(forms.Form):

    def __init__(self, *args, **kwargs):

        analysis = kwargs.pop("analysis")

        super().__init__(*args, **kwargs)

        analysis = util.safe_loads(analysis)
        initial = json.dumps(analysis, indent=4)

        self.fields["text"] = forms.CharField(initial=initial)
        self.fields["save_or_preview"] = forms.CharField(initial="preview")
        # put in function.
        for field in analysis:

            data = analysis[field]
            display_type = data["display_type"]
            factory.check_display(display_type)

            if data.get("visible") == 1:

                self.fields[field] = factory.TYPE2FUNC[display_type](data)

    def save(self, *args, **kwargs):

        super(EditAnalysis, self).save(*args, **kwargs)




# class ResultForm(forms.ModelForm):
#
#     text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))
#
#     class Meta:
#         model = Result
#         fields = ['title', 'text', 'commands', 'state', 'directory']
#
#


