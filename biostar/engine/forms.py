import hjson, logging
from django import forms
from django.utils.safestring import mark_safe
from .models import Project, Data, Analysis, Job, Access
from . import tasks
from .const import *
import os
from . import models, auth


# Share the logger with models.
logger = models.logger

def join(*args):
    return os.path.abspath(os.path.join(*args))


class ProjectForm(forms.ModelForm):
    image = forms.ImageField(required=False)

    class Meta:
        model = Project
        fields = ['name', 'summary', 'text', 'image', "privacy", "sticky"]


class DataUploadForm(forms.ModelForm):
    #choices = DATA_TYPES.items()
    #data_type = forms.IntegerField(widget=forms.Select(choices=choices))

    file = forms.FileField()

    class Meta:
        model = Data
        fields = ['file', 'summary', 'text', "sticky"]


class DataEditForm(forms.ModelForm):
    #choices = DATA_TYPES.items()
    #data_type = forms.IntegerField(widget=forms.Select(choices=choices))

    class Meta:
        model = Data
        fields = ['name', 'summary', 'text', 'sticky']


class AnalysisEditForm(forms.ModelForm):

    class Meta:
        model = Analysis
        fields = ['name', "image",'text', "summary", 'sticky']


class JobEditForm(forms.ModelForm):

    class Meta:
        model = Job
        fields = ['name', "image",'text','sticky']


class GrantAccess(forms.Form):

    users = forms.IntegerField()
    add_or_remove = forms.CharField(initial="")

    def __init__(self, project, current_user, access,*args, **kwargs):
        self.project = project
        self.user = current_user
        self.access = access

        super().__init__(*args, **kwargs)

    def is_valid(self, request=None):
        valid = super(GrantAccess, self).is_valid()

        # Only users with admin privilege or higher get to grant access to projects
        admin_only  = auth.check_obj_access(user=self.user, instance=self.project,
                                            request=request,access=Access.ADMIN_ACCESS)
        return valid and admin_only


    def process(self, add=False,remove=False):

        assert not (add and remove), "Can add or remove, not both"

        # More than one can be selected
        users = self.data.getlist('users')
        added, removed, errmsg  = 0,0, []

        for user_id in users:
            user = models.User.objects.filter(id=user_id).first()
            has_access = user.access_set.filter(project=self.project).first()

            # Can only add people without access
            addcond = (not has_access or has_access.access == Access.PUBLIC_ACCESS)

            # Can only remove people with access
            remcond = (has_access and has_access.access > Access.PUBLIC_ACCESS)

            if add and addcond:

                added += 1
                if not has_access:
                    access = Access.objects.create(user=user, project=self.project, access=self.access)
                    access.save()
                    continue
                has_access.access = self.access
                has_access.save()

            elif remove and remcond:
                # Changes access to Access.PUBLIC_ACCESS
                has_access.access = Access.PUBLIC_ACCESS
                has_access.save()
                removed += 1

            # Trying to add or remove user without meeting conds not allowed
            elif (add and (not addcond)) or (remove and (not remcond)):
                errmsg.append(f"{user.first_name}")

        if errmsg:
            errmsg = f"{', '.join(errmsg)} already in project" if add else \
                f"Can not remove: {', '.join(errmsg)}"

        return added, removed, errmsg


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

            tasks.copier(target_project=self.project.id, fname=path, link=True)

            logger.info(f"Copy data at: {path}")

        return len(paths)


class AnalysisCopyForm(forms.Form):

    projects = forms.IntegerField()

    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis
        super().__init__(*args, **kwargs)

    def process(self):

        projects = self.data.getlist('projects')

        for project_id in projects:
            current_project = Project.objects.filter(id=project_id).first()

            current_params = self.analysis_params(project=current_project)
            new_analysis = auth.create_analysis(**current_params)
            # Images needs to be set by it set
            new_analysis.image.save(self.analysis.name, self.analysis.image, save=True)
            new_analysis.save()

        return len(projects)


    def analysis_params(self, project=None):

        project = project or self.analysis.project
        json_text, template = self.analysis.json_text, self.analysis.template
        owner, summary = self.analysis.owner, self.analysis.summary
        name, text = f"Copy of: {self.analysis.name}", self.analysis.text

        params = dict(project=project, json_text=json_text, template=template,
                      user=owner, summary=summary, name=name, text=text)
        return params


class RunAnalysis(forms.Form):

    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis
        self.json_data = self.analysis.json_data
        self.project = self.analysis.project

        super().__init__(*args, **kwargs)

        # This loop needs to be here to register the fields and trigger is_valid() later on.
        for name, data in self.json_data.items():
            field = auth.make_form_field(data, self.project)
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

            # If it has a path it is an uploaded file.
            if obj.get("path") or obj.get("link"):

                data_id = self.cleaned_data.get(field, '')
                data_id = int(data_id)
                data = datamap.get(data_id)
                data.fill_dict(obj)

            if field in self.cleaned_data:
                obj["value"] = self.cleaned_data[field]
        return json_data



class EditAnalysisForm(forms.Form):

    save_or_preview = forms.CharField(initial="preview")
    def __init__(self, analysis, *args, **kwargs):

        self.analysis = analysis

        super().__init__(*args, **kwargs)

        json_data = hjson.loads(self.analysis.json_text)
        initial = hjson.dumps(json_data, indent=4)

        self.fields["json_text"] = forms.CharField(initial=initial)
        self.fields["template"] = forms.CharField(initial=self.analysis.template)
        self.generate_form(json_data)


    def preview(self):

        json_data = hjson.loads(self.cleaned_data["json_text"].rstrip())

        # Refresh form with most recent json data
        self.generate_form(json_data)


    def save(self):

        super(EditAnalysisForm, self).clean()
        json_data = hjson.loads(self.cleaned_data["json_text"])

        # Refresh form
        self.generate_form(json_data)

        spec = hjson.loads(self.cleaned_data["json_text"])

        if spec.get("settings"):
            self.analysis.name = spec["settings"].get("name", self.analysis.name)
            self.analysis.text = spec["settings"].get("text", self.analysis.text)

        self.analysis.json_text = self.cleaned_data["json_text"]

        #TODO: test more ( probs need to sluggify both)
        if self.analysis.template != self.cleaned_data["template"]:
            self.analysis.security = Analysis.UNDER_REVIEW

        self.analysis.template = self.cleaned_data["template"]

        self.analysis.save()

        return self.analysis


    def generate_form(self, json_obj):

        for name, obj in json_obj.items():
            field = auth.make_form_field(obj, self.analysis.project)

            if field:
                self.fields[name] = field
