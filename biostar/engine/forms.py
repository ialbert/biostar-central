import copy
import shlex
import hjson
import io
import os
from django import forms
from django.db.models import Sum
from django.utils.safestring import mark_safe
from django.contrib import messages
from django.urls import reverse
from django.conf import settings

from biostar.accounts.models import User, Profile
from . import models, auth, factory, util
from .const import *
from .models import Project, Data, Analysis, Job, Access

# Share the logger with models.
logger = models.logger

TEXT_UPLOAD_MAX = 10000


def join(*args):
    return os.path.abspath(os.path.join(*args))


def check_size(fobj, maxsize=0.3):
    # maxsize in megabytes!

    try:
        if fobj and fobj.size > maxsize * 1024 * 1024.0:
            curr_size = fobj.size / 1024 / 1024.0
            msg = f"File too large: {curr_size:0.1f}MB should be < {maxsize:0.1f}MB"
            raise forms.ValidationError(msg)
    except Exception as exc:
        raise forms.ValidationError(f"File size validation error: {exc}")

    return fobj


def check_upload_limit(file, user):
    "Checks if the intended file pushes user over their upload limit."

    uploaded_files = Data.objects.filter(owner=user, method=Data.UPLOAD)
    currect_size = uploaded_files.aggregate(Sum("size"))["size__sum"] or 0
    to_mb = lambda x: x / 1024 / 1024

    projected = to_mb(file.size + currect_size)
    max_mb = user.profile.max_upload_size

    allowed = (max_mb - to_mb(currect_size)) or 0

    if projected > max_mb:
        msg = f"<b>Over your {max_mb:0.0001f} MB total upload limit.</b> "
        msg = msg + f"""
                File too large: currently <b>{to_mb(file.size):0.0001f} MB</b>
                should be <b> < {allowed:0.001f} MB</b>
                """
        raise forms.ValidationError(mark_safe(msg))

    return file


def clean_file(fobj, user, project, check_name=True):

    if not fobj:
        return fobj

    check_size(fobj=fobj, maxsize=settings.MAX_FILE_SIZE_MB)
    check_upload_limit(file=fobj, user=user)

    # Check if this name already exists.
    if check_name and Data.objects.filter(name=fobj.name, project=project).exists():
        msg = "Name already exists. Upload another file or rename existing data."
        raise forms.ValidationError(msg)

    return fobj


class ProjectForm(forms.ModelForm):
    image = forms.ImageField(required=False)

    # Should not edit uid because data directories get recreated
    # uid = forms.CharField(max_length=32, required=False)
    choices = list(filter(lambda x: x[0] != Project.SHAREABLE, Project.PRIVACY_CHOICES))
    privacy = forms.IntegerField(widget=forms.Select(choices=choices))

    class Meta:
        model = Project
        fields = ['name', 'summary', 'text', 'image', "privacy", "sticky"]

    def clean_image(self):
        cleaned_data = super(ProjectForm, self).clean()
        image = cleaned_data.get('image')
        check_size(fobj=image)

        return image

    def custom_save(self, owner):
        """Used to save on creation using custom function."""

        name = self.cleaned_data["name"]
        text = self.cleaned_data["text"]
        summary = self.cleaned_data["summary"]
        stream = self.cleaned_data["image"]
        sticky = self.cleaned_data["sticky"]
        privacy = self.cleaned_data["privacy"]
        project = auth.create_project(user=owner, name=name, summary=summary, text=text,
                                      stream=stream, sticky=sticky, privacy=privacy)
        project.save()

        return project


class DataUploadForm(forms.ModelForm):

    file = forms.FileField(required=False)
    input_text = forms.CharField(max_length=TEXT_UPLOAD_MAX, required=False)
    data_name = forms.CharField(required=False)
    type = forms.CharField(max_length=32, required=False)

    def __init__(self, user, project, *args, **kwargs):
        self.user = user
        self.project = project
        super().__init__(*args, **kwargs)

    def save(self, **kwargs):

        text = self.cleaned_data["text"]
        stream = self.cleaned_data["file"]
        input_text = self.cleaned_data['input_text']
        summary = self.cleaned_data["summary"]
        type = self.cleaned_data["type"]
        name = self.cleaned_data['data_name']

        if stream:
            name = name or stream.name
        else:
            stream = io.StringIO(initial_value=input_text)

        data = auth.create_data(stream=stream, name=name, text=text, user=self.user,
                                project=self.project, summary=summary, type=type)
        if input_text and not self.cleaned_data["file"]:
            Data.objects.filter(pk=data.pk).update(method=Data.TEXTAREA)
            stream.close()

        return data

    class Meta:
        model = Data
        fields = ['data_name', 'file', 'input_text', 'summary', 'text', "sticky", "type"]

    def clean(self):
        cleaned_data = super(DataUploadForm, self).clean()

        if not (cleaned_data.get("file") or cleaned_data.get("input_text")):
            raise forms.ValidationError("Upload a file or write into the text field to create some data.")

        if cleaned_data.get("input_text") and not cleaned_data.get("file"):
            if not cleaned_data.get("data_name"):
                raise forms.ValidationError("Name is required with text inputs.")
        return cleaned_data

    def clean_file(self):
        cleaned_data = super(DataUploadForm, self).clean()

        return clean_file(fobj=cleaned_data.get('file'),
                          user=self.user,
                          project=self.project, check_name=False)

    def clean_type(self):
        cleaned_data = super(DataUploadForm, self).clean()
        fobj = cleaned_data.get('file')
        if fobj:
            name = fobj.name
        else:
            name = cleaned_data.get('data_name')

        root, ext = os.path.splitext(name)
        ext = ext[1:]
        datatype = EXT_TO_TYPE.get(ext, cleaned_data.get('type'))

        datatype = datatype.upper() or ext.upper()

        return datatype


class DataEditForm(forms.ModelForm):

    type = forms.CharField(max_length=32, required=False)

    def __init__(self, user, *args, **kwargs):

        self.user = user

        super().__init__(*args, **kwargs)

        if self.instance.method == Data.UPLOAD:
            self.fields["file"] = forms.FileField(required=False)

        elif self.instance.method == Data.TEXTAREA:
            initial = ''.join(open(self.instance.get_files()[0], 'r').readlines())
            self.fields["input_text"] = forms.CharField(max_length=TEXT_UPLOAD_MAX,
                                                        required=True,
                                                        initial=initial)

    def save(self, commit=True):

        cleaned_data = super(DataEditForm, self).clean()
        fobj = cleaned_data.get('file')
        input_text = cleaned_data.get("input_text")
        current_file = self.instance.get_files()[0]

        if input_text:
            fobj = io.StringIO(initial_value=input_text)

        if fobj:
            util.write_stream(stream=fobj, dest=current_file)

        return super(DataEditForm, self).save(commit)


    class Meta:
        model = Data
        fields = ['name', 'summary', 'text', 'sticky', "type"]


    def clean_file(self):
        cleaned_data = super(DataEditForm, self).clean()
        return clean_file(fobj=cleaned_data.get('file'),
                          user=self.user,
                          project=self.instance.project,
                          check_name=False)

    def clean_type(self):
        cleaned_data = super(DataEditForm, self).clean()

        datatype = cleaned_data.get('type')

        datatype = datatype.upper()

        return datatype


class RecipeCodeEdit(forms.ModelForm):

    uid = forms.CharField(max_length=32, required=False)

    def __init__(self, user, recipe, *args, **kwargs):
        self.user = user
        self.recipe = recipe

        super().__init__(*args, **kwargs)

    class Meta:
        model = Analysis
        fields = ["template"]

    def clean(self):

        # Check if the user is a manager or has write access before making changes.
        entry = auth.check_obj_access(user=self.user, instance=self.recipe, access=Access.WRITE_ACCESS,
                                      login_required=True, role=Profile.MANAGER)
        if not entry:
            raise forms.ValidationError("You need write access to change the code.")

    # Turn all input into Unix line ending.
    def clean_template(self):
        cleaned_data = super(RecipeCodeEdit, self).clean()
        template = cleaned_data.get('template')
        template = "\n".join(template.splitlines())
        return template


class RecipeForm(forms.ModelForm):
    image = forms.ImageField(required=False)
    uid = forms.CharField(max_length=32, required=False)


    class Meta:
        model = Analysis
        fields = ["name", "sticky", "image", "summary", "text", "uid"]

    def clean_image(self):
        cleaned_data = super(RecipeForm, self).clean()
        image = cleaned_data.get('image')
        check_size(fobj=image)

        return image

    def clean_uid(self):
        cleaned_data  =  super(RecipeForm, self).clean()
        uid = cleaned_data.get('uid')
        if uid and not uid.isalnum():
            msg = "Only alphanumeric characters allowed, no spaces."
            raise forms.ValidationError(msg)

        return uid


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
        "Changes users access to a project"

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

        return user, new_access


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


class RecipeDiff(forms.Form):

    REVERT, APPROVE = 'REVERT', 'APPROVE'

    action = forms.CharField()

    def __init__(self, recipe, user, request, *args, **kwargs):

        self.recipe = recipe
        self.user = user
        self.request = request
        super().__init__(*args, **kwargs)

    def save(self):

        action = self.cleaned_data.get("action")

        if action == self.REVERT:
            self.recipe.template = self.recipe.last_valid
            messages.success(self.request, "Recipe has been reverted to original.")

        elif action == self.APPROVE:
            self.recipe.last_valid = self.recipe.template
            messages.success(self.request, "Recipe changes have been approved.")
        self.recipe.security = Analysis.AUTHORIZED
        self.recipe.save()
        return self.recipe

    def clean(self):

        cleaned_data = super(RecipeDiff, self).clean()
        action = cleaned_data.get("action")
        msg = "You don't have sufficient access rights to overwrite this entry."
        has_access = auth.check_obj_access(user=self.user, request=self.request, instance=self.recipe,
                                           role=Profile.MANAGER, access=Access.WRITE_ACCESS, login_required=True)
        if not has_access:
            raise forms.ValidationError(msg)

        if action not in (self.REVERT, self.APPROVE):
            raise forms.ValidationError("Can only revert or approve changes.")

        if action == self.APPROVE and not self.user.profile.is_manager:
            msg = "You have to be a manager."
            raise forms.ValidationError(msg)


class RecipeInterface(forms.Form):
    # The name of results when running the recipe.
    # name = forms.CharField(max_length=256, label="Name", help_text="This is how you can identify the run.")

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

        # The owner of the project may run any recipe in the project.
        if self.project.owner == self.user:
            return

        # Check the permissions for the recipe.
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
