import hjson
from django import forms
import hjson
from . import models, auth
from . import tasks
from .const import *
from .models import Project, Data, Analysis, Job, Access

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
    # choices = DATA_TYPES.items()
    # data_type = forms.IntegerField(widget=forms.Select(choices=choices))

    file = forms.FileField()

    class Meta:
        model = Data
        fields = ['file', 'summary', 'text', "sticky"]


class DataEditForm(forms.ModelForm):
    # choices = DATA_TYPES.items()
    # data_type = forms.IntegerField(widget=forms.Select(choices=choices))

    class Meta:
        model = Data
        fields = ['name', 'summary', 'text', 'sticky']


class AnalysisEditForm(forms.ModelForm):
    class Meta:
        model = Analysis
        fields = ['name', "image", 'text', "summary", 'sticky']


class JobEditForm(forms.ModelForm):
    class Meta:
        model = Job
        fields = ['name', "image", 'text', 'sticky']


class GrantAccess(forms.Form):
    users = forms.IntegerField()

    def __init__(self, project, current_user, access, *args, **kwargs):
        self.project = project
        self.user = current_user
        self.access = access

        super().__init__(*args, **kwargs)

    def is_valid(self, request=None):
        valid = super(GrantAccess, self).is_valid()

        # Only users with admin privilege or higher get to grant access to projects
        admin_only = auth.check_obj_access(user=self.user, instance=self.project,
                                           request=request, access=Access.ADMIN_ACCESS)
        return valid and admin_only

    def process(self, add=False, remove=False):

        assert not (add and remove), "Can add or remove, not both"

        # More than one can be selected
        users = self.data.getlist('users')
        added, removed, errmsg = 0, 0, []

        for user_id in users:
            user = models.User.objects.filter(id=user_id).first()
            has_access = user.access_set.filter(project=self.project).first()

            # Can only add people without access
            addcond = (not has_access or has_access.access == Access.NO_ACCESS)

            # Can only remove people with access
            remcond = (has_access and has_access.access > Access.NO_ACCESS)

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
                has_access.access = Access.NO_ACCESS
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

    # TODO: refractor asap; does not need to be a list only one is picked
    def process(self):
        projects = self.data.getlist('projects')
        project_id = projects[0]

        if project_id == "0":
            return projects, None

        current_project = Project.objects.filter(id=project_id).first()

        current_params = auth.get_analysis_attr(analysis=self.analysis, project=current_project)
        new_analysis = auth.create_analysis(**current_params)

        # Images needs to be set by it set
        new_analysis.image.save(self.analysis.name, self.analysis.image, save=True)
        new_analysis.name = f"Copy of: {self.analysis.name}"
        new_analysis.state = self.analysis.state
        new_analysis.security = self.analysis.security
        new_analysis.save()

        return projects, new_analysis


class NameInput(forms.TextInput):
    template_name = 'interface/name.html'


class RecipeInterface(forms.Form):
    name = forms.CharField(max_length=256, widget=NameInput)

    def __init__(self, project, json_data, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # The json data determines what fields does the form have.
        self.json_data = json_data

        # The project is required to select data from.
        self.project = project

        # Insert the dynamic fields into the form.
        for name, data in self.json_data.items():
            field = auth.make_form_field(data, self.project)
            if field:
                self.fields[name] = field

    def process(self):
        '''
        Replaces the value of data fields with the path to the data.
        Should be called after the form has been filled and is valid.
        '''
        # Gets all data for the project
        datamap = dict((data.id, data) for data in self.project.data_set.all())

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


class EditCode(forms.Form):
    SAVE = "SAVE"

    # Determines what action to perform on the form.
    action = forms.CharField()

    # The script template.
    template = forms.CharField()

    # The json specification.
    json = forms.CharField()

    def __init__(self, user, project, *args, **kwargs):
        self.user = user
        self.project = project
        super().__init__(*args, **kwargs)

    def clean_json(self):
        cleaned_data = super(EditCode, self).clean()
        json_text = cleaned_data.get("json")
        try:
            hjson.loads(json_text)
        except Exception as exc:
            msg = f"Invalid json: {exc}"
            raise forms.ValidationError(msg)
        return json_text

    def clean(self):
        cleaned_data = super(EditCode, self).clean()
        action = cleaned_data.get("action")

        if action == self.SAVE:
            msg = "You don't have sufficient access rights to overwrite this entry."
            if self.user.is_anonymous():
                raise forms.ValidationError(msg)
            entry = Access.objects.filter(user=self.user, project=self.project).first()
            if not entry or entry.access < Access.EDIT_ACCESS:
                raise forms.ValidationError(msg)
