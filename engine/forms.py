import hjson as json

from django import forms
from django.utils.translation import gettext_lazy as helpers
from django.contrib.auth.models import User
from .models import Project, Data, Analysis
from . import util
from . import factory

from pagedown.widgets import PagedownWidget
from engine.const import *
from engine.web.auth import get_data

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

    #text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))

    class Meta:

        model = Project
        fields = ['title', 'text']


class DataForm(forms.Form):

    text = forms.CharField(widget=forms.Textarea(), max_length=512)
    file = forms.FileField(label="Upload Data File")
    type = forms.IntegerField(widget=forms.Select(choices=Data.TYPE_CHOICES),
                              initial=Data.FILE)

    def save(self, *args, **kwargs):

        super(DataForm, self).save(*args, **kwargs)


def make_form_field(value):

    visible = value.get("visible")
    display_type = value.get("display_type", '')

    if not display_type:
        return ""

    factory.check_display(display_type)

    if visible == 1:
        return factory.TYPE2FUNC[display_type](value)

    return  ""


class RunAnalysis(forms.Form):

    def __init__(self, *args, **kwargs):

        self.json_data = kwargs.pop("json_data")
        self.json_data = json.loads(self.json_data)

        super().__init__(*args, **kwargs)

        # Job needs a title
        self.fields["title"] = forms.CharField(max_length=256,
                                               help_text="Leave Empty to fill with Analysis name.",
                                               required=False)

        for field, value in self.json_data.items():
            form_field = make_form_field(value)

            if form_field :
                self.fields[field] = form_field


    def save(self, *args, **kwargs):

        super(RunAnalysis, self).save(*args, **kwargs)


    def process(self):
        '''
        Replaces the value of data fields with the path to the data.
        Should be called after the form has been filled and is valid.
        '''

        # Should make a copy of the json
        for field, obj in self.json_data.items():
            # No need to read fields that were not set.
            if not obj.get(FIELD_VISIBLE):
                continue

            if obj.get(FIELD_ORIGIN) == PROJECT_ORIGIN:
                data_id = self.cleaned_data.get(field, 0)
                data = get_data(data_id)
                obj["value"] = data.get_path().path
                continue
            if field in self.cleaned_data:
                # Mutates the value key.
                obj["value"] = self.cleaned_data[field]


        print (json.dumps(self.json_data))

        return self.json_data


class EditAnalysis(forms.Form):

    def __init__(self, *args, **kwargs):

        analysis = kwargs.pop("analysis")

        super().__init__(*args, **kwargs)

        analysis = util.safe_loads(analysis)
        initial = json.dumps(analysis, indent=4)

        self.fields["text"] = forms.CharField(initial=initial)
        self.fields["save_or_preview"] = forms.CharField(initial="preview")


        for field, value in analysis.items():
            form_field = make_form_field(value)

            if form_field:
                self.fields[field] = form_field

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


