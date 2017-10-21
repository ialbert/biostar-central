import hjson as json

from django import forms
from django.utils.translation import gettext_lazy as helpers
from django.contrib.auth.models import User, Group
from .models import Project, Data, Analysis
from . import util
from . import factory
from pagedown.widgets import PagedownWidget
from engine.const import *
from engine.web.auth import get_data
from . import models

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
    class Meta:
        model = Project
        fields = ['title', 'summary', 'text']


class DataUploadForm(forms.ModelForm):
    class Meta:
        model = Data
        fields = ['file', 'summary', 'text']


class DataEditForm(forms.ModelForm):
    class Meta:
        model = Data
        fields = ['name', 'summary', 'text']


def make_field(obj, project):
    field = ''
    visible = obj.get("visible")
    path = obj.get("path")
    display_type = obj.get("display_type", '')

    if not display_type:
        return field

    if visible:
        if path:
            data_type = obj.get("data_type")
            field = factory.data_generator(obj, project=project, data_type=data_type)
        else:
            field = factory.TYPE2FUNC[display_type](obj)

    return field



class ExportAnalysis(forms.Form):


    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis
        # Was told to include all projects
        #projects = [(proj.id, proj.title) for proj in Project.objects.all() if proj.id!=self.analysis.project.id]
        projects = [(proj.id, proj.title) for proj in Project.objects.all()]
        super().__init__(*args, **kwargs)
        self.fields["project"] = forms.IntegerField(widget=forms.Select(choices=projects))


    def export(self):
        exported_to = self.cleaned_data.get("project")
        project = Project.objects.filter(id=exported_to).first()

        json_text, template = self.analysis.json_text, self.analysis.template
        owner, summary = self.analysis.owner, self.analysis.summary
        title, text = self.analysis.title, self.analysis.text

        analysis = project.create_analysis(json_text, template, owner, summary, title, text)
        analysis.save()
        return project, analysis


class RunAnalysis(forms.Form):
    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis
        self.json_data = json.loads(self.analysis.json_text)

        super().__init__(*args, **kwargs)

        self.fields["title"] = forms.CharField(max_length=256,
                                               help_text="Results Title",
                                               required=False)

        for name, obj in self.json_data.items():
            field = make_field(obj, analysis.project)
            # print(name, obj)
            if field:
                self.fields[name] = field

    def save(self, *args, **kwargs):
        super(RunAnalysis, self).save(*args, **kwargs)

    def process(self):
        '''
        Replaces the value of data fields with the path to the data.
        Should be called after the form has been filled and is valid.
        '''
        project = self.analysis.project

        # Gets all data for the project
        datamap = project.get_data()

        print ( type(list(datamap.keys())[0]) )

        json_data = self.json_data.copy()

        for field, obj in json_data.items():
            # No need to read fields that were not set.
            if not obj.get(FIELD_VISIBLE):
                continue

            # If it has a path it is an uploaded file.
            if obj.get("path"):
                data_id = self.cleaned_data.get(field, '')
                data_id = int(data_id)
                data = datamap.get(data_id)
                data.fill_dict(obj)
                continue

            if field in self.cleaned_data:
                # Mutates the value key.
                obj["value"] = self.cleaned_data[field]

            # TODO CHANGE
            if obj.get("path"):
                obj["path"] = os.path.abspath(obj.get("path"))

        return json_data


class EditAnalysisForm(forms.Form):
    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis

        super().__init__(*args, **kwargs)

        json_data = json.loads(self.analysis.json_text)
        initial = json.dumps(json_data, indent=4)

        self.fields["text"] = forms.CharField(initial=initial)
        self.fields["save_or_preview"] = forms.CharField(initial="preview")

        self.generate_form(json_data)

    def preview(self):

        cleaned_data = super(EditAnalysisForm, self).clean()
        json_data = json.loads(cleaned_data["text"])

        self.generate_form(json_data)

    def save(self):

        cleaned_data = super(EditAnalysisForm, self).clean()
        json_data = json.loads(cleaned_data["text"])

        self.generate_form(json_data)

        spec = json.loads(self.cleaned_data["text"])
        filler = dict(display_type='')

        if spec.get("settings", filler).get("display_type") == "MODEL":
            self.analysis.title = spec["settings"].get("title", self.analysis.title)
            self.analysis.text = spec["settings"].get("text", self.analysis.text)

        self.analysis.json_text = self.cleaned_data["text"]
        self.analysis.save()

        return self.analysis

    def save_to_file(self, file):
        return

    def generate_form(self, json_obj):

        for name, obj in json_obj.items():
            field = make_field(obj, self.analysis.project)

            if field:
                self.fields[name] = field
