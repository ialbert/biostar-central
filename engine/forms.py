import json

from django import forms
from django.utils.translation import gettext_lazy as helpers
from django.contrib.auth.models import User
from .models import Project, Data, Analysis
from .util import safe_load, TYPE2FUNC, handle_no_type

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


class DataForm(forms.ModelForm):

    text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))
    class Meta:
        model = Data
        fields = ['title', 'text']

    def cleaned_data(self, *args):
        return self.cleaned_data


class AnalysisForm(forms.ModelForm):

    json_spec = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))
    class Meta:
        model = Analysis
        fields = ['title', "json_spec"]


class RunForm(forms.Form):


    def __init__(self, *args, **kwargs):

        json_spec = kwargs.pop("json_spec")

        super().__init__(*args, **kwargs)

        json_spec = safe_load(json_spec)

        for field in json_spec:

            field_template = TYPE2FUNC.get(field["type"], handle_no_type)(field)

            if field.get("visible") == 1:
                exec(field_template)


    def save(self, *args, **kwargs):

        json_spec = {}
        for f in self.fields:
            if "json_" in f:
                json_spec[f] = self.fields[f]

        self.fields["json_spec"] = json.dumps(json_spec)

        super(RunForm, self).save(*args, **kwargs)


class EditForm(forms.Form):


    def __init__(self, *args, **kwargs):

        json_spec = kwargs.pop("json_spec")

        super().__init__(*args, **kwargs)

        self.json_spec = json_spec

        json_spec = safe_load(json_spec)
        self.fields["text"] = forms.CharField(initial=json.dumps(json_spec, indent=4))

        for field in json_spec:

            field_template = TYPE2FUNC.get(field["type"], handle_no_type)(field)

            if field.get("visible") == 1:
                exec(field_template)


    def save(self, *args, **kwargs):

        json_spec = {}
        # write the text feild into the json_spec file ( thats the update).
        for f in self.fields:
            if "json_" in f:
                json_spec[f] = self.fields[f]

        self.fields["json_spec"] = json.dumps(json_spec)

        super(EditForm, self).save(*args, **kwargs)


    def preview(self):
        # change the
        return









# class ResultForm(forms.ModelForm):
#
#     text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))
#
#     class Meta:
#         model = Result
#         fields = ['title', 'text', 'commands', 'state', 'directory']
#
#


