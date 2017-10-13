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


class LogoutForm(forms.Form):
    pass


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


def make_field(obj, project):

    field = ''
    visible = obj.get("visible")
    origin = obj.get(FIELD_ORIGIN)
    display_type = obj.get("display_type", '')

    if not display_type:
        return field

    if visible:
        if origin == PROJECT_ORIGIN:
            data_type = obj.get("data_type")
            field = factory.data_generator(obj, project=project, data_type=data_type)
        else:
            field = factory.TYPE2FUNC[display_type](obj)

    return field


class RunAnalysis(forms.Form):

    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis
        self.json_data = json.loads(self.analysis.json_data)

        super().__init__(*args, **kwargs)

        # Job needs a title
        self.fields["title"] = forms.CharField(max_length=256,
                                               help_text="Leave Empty to fill with Analysis name.",
                                               required=False)

        for name, obj in self.json_data.items():
            field = make_field(obj, analysis.project)

            if field:
                self.fields[name] = field


    def save(self, *args, **kwargs):

        super(RunAnalysis, self).save(*args, **kwargs)


    def process(self):
        '''
        Replaces the value of data fields with the path to the data.
        Should be called after the form has been filled and is valid.
        '''

        # TODO: should make a copy of the json rather than mutated it.
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

        return self.json_data


class EditAnalysis(forms.Form):

    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis
        self.json_data =  kwargs.get('json_data') or self.analysis.json_data
        self.json_data = util.safe_loads(self.json_data)

        if kwargs.get("json_data"):
            kwargs.pop('json_data')

        super().__init__(*args, **kwargs)

        initial = json.dumps(self.json_data, indent=4)

        self.fields["text"] = forms.CharField(initial=initial)
        self.fields["save_or_preview"] = forms.CharField(initial="preview")

        for name, obj in self.json_data.items():
            field = make_field(obj, analysis.project)

            if field:
                self.fields[name] = field

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


