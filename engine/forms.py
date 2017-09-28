import json

from django import forms
from django.utils.translation import gettext_lazy as helpers
from django.contrib.auth.models import User
from .models import Project, Data, Analysis, Result
from django.forms.extras.widgets import SelectDateWidget
from .util import *

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

    text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))
    class Meta:
        model = Analysis
        fields = ['title', 'text']


class RunForm(forms.Form):

    def __init__(self, *args, **kwargs):

        analysis_id = kwargs.pop("id")
        json_spec = kwargs.pop("json_spec")
        self.makefile = kwargs.pop("makefile")

        super().__init__(*args, **kwargs)

        analysis = Analysis.objects.filter(id=analysis_id).first()
        DATA = analysis.project.data_set.all()

        CHOICES = []

        for i,d in enumerate(DATA):
            CHOICES.append((i, d))

        self.fields["data"] = forms.CharField(label="Sequencing Data",
                                              widget=forms.Select(choices=CHOICES))
        # json_spec is already a charfeild


        json_spec = json.loads(json_spec)

        for f in json_spec:
            print(f)
            widget = f["widget"]
            form_type = f["form_type"]
            label = f["label"]
            choices = list(f["choices"].items())

            # Name is given json_ prefix to make repopulating json_specs easier.
            exec(f"""
                self.fields['json_'+'{f['name']}'] = forms.{form_type}(
                                                            widget=foms.{widget}(choices={choices}),
                                                            label={label})                                                                                                                       
                 """)

    def save(self, *args, **kwargs):
        # self.fields["json_spec"] = json_spec

        json_spec = {}
        for f in self.fields:
            if "json_" in f:
                json_spec[f] = self.fields[f]

        self.fields["json_spec"] = repr(json_spec)
        self.fields["makefile"] = fill_makefile(repr(json_spec), self.makefile)

        super(RunForm, self).save(*args, **kwargs)

# class ResultForm(forms.ModelForm):
#
#     text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))
#
#     class Meta:
#         model = Result
#         fields = ['title', 'text', 'commands', 'state', 'directory']
#
#


