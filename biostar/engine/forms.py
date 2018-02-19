import copy
import shlex
import hjson

from django import forms
from django.db.models import Sum
from biostar import  settings
from biostar.accounts.models import User
from django.utils.safestring import mark_safe

from . import models, auth, factory
from .const import *
from .models import Project, Data, Analysis, Job, Access


# Share the logger with models.
logger = models.logger


def join(*args):
    return os.path.abspath(os.path.join(*args))


def check_size(fobj, maxsize=0.3):
    # maxsize in megabytes!

    try:

        if fobj and fobj.size > maxsize * 1024 * 1024.0 :
            curr_size = fobj.size / 1024 / 1024.0
            msg = f"File too large: {curr_size:0.1f}MB should be < {maxsize:0.1f}MB"
            raise forms.ValidationError(msg)
    except Exception as exc:
        raise forms.ValidationError(f"File size validation error: {exc}")

    return fobj


def check_upload_limit(file, user):
    "Checks if the intended file pushes user over their upload limit."

    uploaded_files = Data.objects.filter(owner=user, method=Data.UPLOAD)
    currect_size = uploaded_files.aggregate(Sum("size"))

    projected = file.size + currect_size["size__sum"]
    max_mb = settings.MAX_AGG_UPLOAD / 1024 / 1024
    curr_size = file.size / 1024 / 1024
    avail = lambda x: 0 if x < 0 else x

    if projected > settings.MAX_AGG_UPLOAD:

        allowed = avail(max_mb-currect_size["size__sum"])
        msg = f"<b>Over the {max_mb:0.001f} MB total upload limit.</b> "
        msg = msg + f"""
                Your have already uploaded {currect_size["size__sum"]/ 1024/ 1024:0.001f} MB. 
                File too large: currently {curr_size:0.0001f} MB
                should be < {allowed:0.001f} MB
                """
        raise forms.ValidationError(mark_safe(msg))

    return file


class ProjectForm(forms.ModelForm):

    image = forms.ImageField(required=False)

    # Should not edit uid because data directories get recreated
    #uid = forms.CharField(max_length=32, required=False)
    choices = list(filter(lambda x: x[0]!= Project.SHAREABLE, Project.PRIVACY_CHOICES))
    privacy = forms.IntegerField(widget=forms.Select(choices=choices))

    class Meta:
        model = Project
        fields = ['name', 'summary', 'text','image', "privacy", "sticky"]

    def clean_image(self):
        cleaned_data = super(ProjectForm, self).clean()
        image = cleaned_data.get('image')
        check_size(fobj=image)

        return image


class DataUploadForm(forms.ModelForm):

    def __init__(self, user, *args, **kwargs):
        self.user = user
        super().__init__(*args, **kwargs)


    file = forms.FileField(required=True)
    type = forms.CharField(max_length=32, required=False)

    class Meta:
        model = Data
        fields = ['file', 'summary', 'text', "sticky",  "type"]


    def clean_file(self):
        cleaned_data = super(DataUploadForm, self).clean()
        fobj = cleaned_data.get('file')
        check_size(fobj=fobj, maxsize=25)
        #check_upload_limit(file=fobj, user=self.user)
        return fobj


class DataEditForm(forms.ModelForm):
    type = forms.CharField(max_length=32, required=False)

    class Meta:
        model = Data
        fields = ['name', 'summary', 'text', 'sticky', "type"]


class RecipeForm(forms.ModelForm):
    image = forms.ImageField(required=False)

    class Meta:
        model = Analysis
        fields = ["name", "sticky", "image", "summary", "text" ]

    def clean_image(self):
        cleaned_data = super(RecipeForm, self).clean()
        image = cleaned_data.get('image')
        check_size(fobj=image)

        return image


class JobEditForm(forms.ModelForm):

    class Meta:
        model = Job
        fields = ['name', "image", 'text', 'sticky']


class ChangeUserAccess(forms.Form):


    user_id = forms.IntegerField(required=True, widget=forms.HiddenInput())
    project_uid = forms.CharField(required=True, widget=forms.HiddenInput())
    choices = filter(lambda x: x[0] != Access.OWNER_ACCESS, Access.ACCESS_CHOICES)
    access = forms.IntegerField(initial=Access.NO_ACCESS,
                                widget=forms.Select(choices=choices))


    def save(self):
        "Change users access to a project"

        user_id = self.cleaned_data["user_id"]
        project_uid = self.cleaned_data["project_uid"]
        user = User.objects.filter(id=user_id).first()
        project = Project.objects.filter(uid=project_uid).first()
        current = Access.objects.filter(user=user, project=project)

        if current:
            current.update(access=self.cleaned_data.get("access", current.first().access))
            return user, current.first()
        new_access = Access(user=user, project=project,
                            access=self.cleaned_data.get("access", Access.NO_ACCESS))
        new_access.save()

        return user,new_access


def access_forms(users, project, exclude=()):
    """Generate a list of forms for a given user list
    Param exclude: a list of users to exclude"""

    forms = []
    for user in users:
        if user in exclude:
            continue
        access = Access.objects.filter(user=user, project=project).first()
        initial = dict(access=Access.NO_ACCESS, user_id=user.id)
        if access:
            initial = dict(access=access.access, user_id=user.id)

        access_form = ChangeUserAccess(initial=initial)
        forms.append((user, access_form))

    return forms


def clean_text(textbox):
    return shlex.quote(textbox)


class RecipeInterface(forms.Form):

    # The name of results when running the recipe.
    #name = forms.CharField(max_length=256, label="Name", help_text="This is how you can identify the run.")

    def __init__(self, request, analysis, json_data, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # The json data determines what fields does the form have.
        self.json_data = json_data

        # The project is required to select data from.
        self.analysis = analysis
        self.project = analysis.project

        # Get request specific information
        self.request = request
        self.user = self.request.user

        # Create the dynamic field from each key in the data.
        for name, data in self.json_data.items():
            field = factory.dynamic_field(data, self.project)

            # Insert only valid fields.
            if field:
                self.fields[name] = field

    def clean(self):
        cleaned_data = super(RecipeInterface, self).clean()

        if self.user.is_anonymous:
            msg1 = "Only logged in users may execute recipes."
            raise forms.ValidationError(msg1)

        # Check the permissions for
        entry = Access.objects.filter(user=self.user, project=self.project).first()

        if not entry or entry.access < Access.WRITE_ACCESS:
            msg2 = "You don't have write access to this project. Copy this analysis to another project."
            raise forms.ValidationError(msg2)

        if self.analysis.security != Analysis.AUTHORIZED:
            msg3 = "The recipe has been modified. It must be reviewed by a staff member to authorize it."
            raise forms.ValidationError(msg3)

    def fill_json_data(self):
        """
        Produces a filled in JSON data based on user input.
        Should be called after the form has been filled and is valid.
        """

        # Creates a data.id to data mapping.
        store = dict((data.id, data) for data in self.project.data_set.all())

        # Make a copy of the original json data used to render the form.
        json_data = copy.deepcopy(self.json_data)

        # Alter the json data and fill in the extra information.
        for field, item in json_data.items():

            # If the field is a data field then fill in more information.
            if item.get("source") == "PROJECT":
                data_id = int(self.cleaned_data.get(field))
                data = store.get(data_id)
                # This mutates the `item` dictionary!
                data.fill_dict(item)
                continue

            # The JSON value will be overwritten with the selected field value.
            if field in self.cleaned_data:
                value = self.cleaned_data[field]
                # Clean the textbox value
                item["value"] = value if item['display'] != TEXTBOX else clean_text(value)

        return json_data


class DataCopyForm(forms.Form):
    "Used to store multiple data uids in data_clipboard with checkbox inputs"

    def __init__(self, request, *args, **kwargs):
        self.request = request
        super().__init__(*args, **kwargs)

    uids = forms.CharField(max_length=models.MAX_TEXT_LEN)

    def save(self):
        for uid in self.cleaned_data:
            auth.load_data_clipboard(uid=uid, request=self.request)
        return len(self.cleaned_data)

    def clean(self):

        uids = self.data.getlist('uids')
        # Override cleaned_data to later access in save()
        self.cleaned_data = []
        for uid in uids:
            if Data.objects.filter(uid=uid).exists:
                self.cleaned_data.append(uid)


class FileCopyForm(forms.Form):
    "Used to save paths found in jobs/data into files_clipboard"

    def __init__(self, request, uid, root_dir, *args, **kwargs):
        self.uid = uid
        self.request = request
        self.root_dir = root_dir
        super().__init__(*args, **kwargs)

    paths = forms.CharField(max_length=models.MAX_FIELD_LEN)

    def save(self):
        # Save the selected files in clipboard,
        # Note: instance.uid is appended and later used to validate where copied files came from
        self.cleaned_data.append(self.uid)
        self.request.session["files_clipboard"] = self.cleaned_data
        return len(self.cleaned_data[:-1])

    def clean(self):
        paths = self.data.getlist('paths')
        for p in paths:
            if not os.path.exists(join(self.root_dir, p)):
                raise forms.ValidationError(f"{p} does not exist")

        # Override cleaned_data to later access in save()
        self.cleaned_data = list(map(join, [self.root_dir]*len(paths), paths))


class EditCode(forms.Form):
    SAVE = "SAVE"

    # Determines what action to perform on the form.
    action = forms.CharField()

    # The script template.
    template = forms.CharField(required=False)

    # The json specification.
    json = forms.CharField(required=False)

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
            if self.user.is_anonymous:
                raise forms.ValidationError(msg)
            entry = Access.objects.filter(user=self.user, project=self.project).first()
            if not entry or entry.access < Access.WRITE_ACCESS:
                raise forms.ValidationError(msg)
