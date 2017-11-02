import hjson, logging

from django import forms
from django.utils.translation import gettext_lazy as helpers
from django.contrib.auth.models import User, Group
from django.core import management
from .models import Project, Data
from .const import *
import os
from . import util
from . import factory
from . import models

# Share the logger with models.
logger = models.logger

def join(*args):
    return os.path.abspath(os.path.join(*args))


class ProjectForm(forms.ModelForm):
    class Meta:
        model = Project
        fields = ['name', 'summary', 'text']


class DataUploadForm(forms.ModelForm):

    choices = DATA_TYPES.items()
    data_type = forms.IntegerField(widget=forms.Select(choices=choices))

    class Meta:
        model = Data
        fields = ['file', 'summary', 'data_type', 'text']


class DataEditForm(forms.ModelForm):
    choices = DATA_TYPES.items()
    data_type = forms.IntegerField(widget=forms.Select(choices=choices))

    class Meta:
        model = Data
        fields = ['name', 'summary', 'data_type','text']


def make_form_field(data, project):

    display_type = data.get("display_type", '')

    # Fields with no display type are not visible.
    if not display_type:
        return

    # Uploaded data is accessed via paths.
    path = data.get("path")

    if path:
        # Project specific data needs a special field.
        data_type = data.get("data_type")
        field = factory.data_generator(data, project=project, data_type=data_type)
    else:
        func = factory.TYPE2FUNC.get(display_type)
        if not func:
            logger.error(f"Invalid display_type={display_type}")
            return
        field = func(data)

    return field


class DataCopyForm(forms.Form):

    paths = forms.CharField(max_length=256)

    def __init__(self, project, job=None, *args, **kwargs):
        self.project = project
        self.job = job
        super().__init__(*args, **kwargs)

    def process(self):
        # More than one can be selected
        paths = self.data.getlist('paths')
        basedir = '' if not self.job else self.job.path

        for path in paths:
            # Figure out the full path based on existing data
            if path.startswith("/"):
                path = path[1:]
            path = join(basedir, path)

            management.call_command("admintasks", copy=True, fname=path, pid=self.project.id)

            print(f"Copy data at: {path}")

        return len(paths)


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

        datamap = dict((data.id, data) for data in self.project.data_set.all() )
        json_data = self.json_data.copy()

        for field, obj in json_data.items():
            # No need to read fields that were not set.

            # If it has a path it is an uploaded file.
            if obj.get("path"):

                data_id = self.cleaned_data.get(field, '')
                data_id = int(data_id)
                data = datamap.get(data_id)
                data.fill_dict(obj)

            if field in self.cleaned_data:
                obj["value"] = self.cleaned_data[field]
        return json_data



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
