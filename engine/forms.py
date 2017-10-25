import hjson, logging

from django import forms
from django.utils.translation import gettext_lazy as helpers
from django.contrib.auth.models import User, Group
from .models import Project, Data, Analysis
from . import util
from . import factory
#from pagedown.widgets import PagedownWidget
from engine.const import *
from engine.web.auth import get_data
from . import models

# Share the logger with models.
logger = models.logger

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
        fields = ['name', 'summary', 'text']


class DataUploadForm(forms.ModelForm):
    class Meta:
        model = Data
        fields = ['file', 'summary', 'text']


class DataEditForm(forms.ModelForm):
    choices = [(y, x) for x,y in DATA_TYPES.items()]
    data_type = forms.IntegerField(widget=forms.Select(choices=choices))
    class Meta:
        model = Data
        fields = ['name', 'summary', 'data_type','text']



def make_form_field(data, project):

    # Should this field be rendered?
    display_type = data.get("display_type", '')
    if not display_type:
        return

    # Is this an existing data
    path = data.get("path")
    origin = data.get('origin')

    # Project specific data
    if path and origin == PROJECT_ORIGIN:
        data_type = data.get("data_type")
        field = factory.data_generator(data, project=project, data_type=data_type)

    # not project specific data
    elif origin == ALL_ORIGIN:
        data_type = data.get("data_type")
        field = factory.data_generator(data, data_type=data_type)

    else:
        func = factory.TYPE2FUNC.get(display_type)
        if not func:
            logger.error(f"Invalid display_type={display_type}")
            return
        field = func(data)

    return field


class ExportData(forms.Form):


    def __init__(self, project, *args, **kwargs):

        self.project = project
        super().__init__(*args, **kwargs)

        self.fields["filename"] = forms.CharField(initial="filename")

    def export(self):

        to_export = self.cleaned_data.get("filename")
        data = self.project.create_data(fname=to_export)
        return data


class ExportAnalysis(forms.Form):


    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis
        projects = [(proj.id, proj.name) for proj in Project.objects.all()]
        super().__init__(*args, **kwargs)
        self.fields["project"] = forms.IntegerField(widget=forms.Select(choices=projects))


    def export(self):
        exported_to = self.cleaned_data.get("project")
        project = Project.objects.filter(id=exported_to).first()

        json_text, template = self.analysis.json_text, self.analysis.template
        owner, summary = self.analysis.owner, self.analysis.summary
        name, text = self.analysis.name, self.analysis.text

        analysis = project.create_analysis(json_text, template, owner, summary, name=name, text=text)
        analysis.save()
        return project, analysis



class RunAnalysis(forms.Form):

    name = forms.CharField(max_length=256, help_text="Results Title", required=True)

    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis
        self.json_data = self.analysis.json_data
        self.project = self.analysis.project

        super().__init__(*args, **kwargs)

        for name, data in self.json_data.items():
            field = make_form_field(data, self.project)
            if field:
                self.fields[name] = field

    def save(self, *args, **kwargs):
        super(RunAnalysis, self).save(*args, **kwargs)

    def process(self):
        '''
        Replaces the value of data fields with the path to the data.
        Should be called after the form has been filled and is valid.
        '''

        # Gets all data for the project

        datamap = self.project.get_data()
        json_data = self.json_data.copy()

        for field, obj in json_data.items():
            # No need to read fields that were not set.
            if not obj.get(FIELD_VISIBLE):
                continue

            # If it has a path it is an uploaded file.
            if obj.get("path") and obj.get("origin")== PROJECT_ORIGIN:

                self.fill_data(field=field, datamap=datamap, obj=obj)
                continue

            elif obj.get("path") and obj.get("origin")== ALL_ORIGIN:

                datamap = dict((obj.id, obj) for obj in Data.objects.all())
                self.fill_data(field=field, datamap=datamap, obj=obj)
                continue

            if field in self.cleaned_data:
                obj["value"] = self.cleaned_data[field]
        return json_data


    def fill_data(self, field, datamap, obj):

        data_id = self.cleaned_data.get(field, '')
        data_id = int(data_id)
        data = datamap.get(data_id)
        data.fill_dict(obj)
        return


class EditAnalysisForm(forms.Form):
    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis

        super().__init__(*args, **kwargs)

        json_data = hjson.loads(self.analysis.json_text)
        initial = hjson.dumps(json_data, indent=4)

        self.fields["text"] = forms.CharField(initial=initial)
        self.fields["save_or_preview"] = forms.CharField(initial="preview")

        self.generate_form(json_data)

    def preview(self):

        cleaned_data = super(EditAnalysisForm, self).clean()
        json_data = hjson.loads(cleaned_data["text"])

        self.generate_form(json_data)

    def save(self):

        cleaned_data = super(EditAnalysisForm, self).clean()
        json_data = hjson.loads(cleaned_data["text"])

        self.generate_form(json_data)

        spec = hjson.loads(self.cleaned_data["text"])
        filler = dict(display_type='')

        if spec.get("settings", filler).get("display_type") == "MODEL":
            self.analysis.name = spec["settings"].get("name", self.analysis.name)
            self.analysis.text = spec["settings"].get("text", self.analysis.text)

        self.analysis.json_text = self.cleaned_data["text"]
        self.analysis.save()

        return self.analysis


    def generate_form(self, json_obj):

        for name, obj in json_obj.items():
            field = make_form_field(obj, self.analysis.project)

            if field:
                self.fields[name] = field
