import json
from string import Template

from django import forms
from django.utils.translation import gettext_lazy as helpers
from django.contrib.auth.models import User
from .models import Project, Data, Analysis

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
        #self.makefile = kwargs.pop("makefile")
        super().__init__(*args, **kwargs)
        json_spec = json.loads(json_spec)

        for f in json_spec:

            choices = f.get("choices")
            min_value, max_value = f.get("min_value"), f.get("max_value")

            # View data already in database
            if f["name"] == "data" and f.get("origin") == "PROJECT":
                # Just loads all data for now ( not project specific ).
                data = Data.objects.all()
                choices = []
                for d in data:
                    choices.append((d.id, d.title))

            # Template used to make that every field.

            template = r"""self.fields['json_'+ '$name']= forms.$form_type(
                                                         widget=forms.$widget(choices=$choices),
                                                         initial='$value',
                                                         label='$label',
                                                         min_value=$min_value,
                                                         max_value=$max_value )""".strip()
            if f.get("visible") == 1:

                data = {"name": f["name"], "form_type": f["form_type"],
                        "label": f["label"], "value": f["value"],
                        "widget": f["widget"]}
                # Some widgets do not have choices
                if choices:
                    data.update(dict(choices=choices))
                else:
                    template = template.replace("choices=$choices", "")

                if min_value and max_value:
                    data.update(dict(min_value=min_value, max_value=max_value))
                else:
                    template = template.replace("min_value=$min_value,", "")
                    template = template.replace("max_value=$max_value", "").strip()

                field_template = Template(template).safe_substitute(data)

                print(field_template)
                #1/0
                exec(field_template)


    def save(self, *args, **kwargs):

        json_spec = {}
        for f in self.fields:
            if "json_" in f:
                json_spec[f] = self.fields[f]

        self.fields["json_spec"] = json.dumps(json_spec)

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


