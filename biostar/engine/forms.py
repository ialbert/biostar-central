import copy
import shlex

import hjson
from django import forms
from django.db.models import Sum
from django.utils.safestring import mark_safe
from django.contrib import messages

from biostar.accounts.models import User, Profile
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

    if projected > max_mb:
        msg = f"<b>Over your {max_mb:0.0001f} MB total upload limit.</b> "
        msg = msg + f"""
                File too large: currently <b>{to_mb(file.size):0.0001f} MB</b>
                should be <b> < {(max_mb-to_mb(currect_size)):0.001f} MB</b>
                """
        raise forms.ValidationError(mark_safe(msg))

    return file


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



class DataUploadForm(forms.ModelForm):
    def __init__(self, user, *args, **kwargs):
        self.user = user
        super().__init__(*args, **kwargs)

    file = forms.FileField(required=True)
    type = forms.CharField(max_length=32, required=False)

    class Meta:
        model = Data
        fields = ['file', 'summary', 'text', "sticky", "type"]

    def clean_file(self):
        cleaned_data = super(DataUploadForm, self).clean()
        fobj = cleaned_data.get('file')

        check_size(fobj=fobj, maxsize=25)
        check_upload_limit(file=fobj, user=self.user)

        return fobj


class DataEditForm(forms.ModelForm):
    type = forms.CharField(max_length=32, required=False)

    class Meta:
        model = Data
        fields = ['name', 'summary', 'text', 'sticky', "type"]


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

        entry = Access.objects.filter(user=self.user, project=self.recipe.project).first()
        entry = entry or Access(access=Access.NO_ACCESS)
        no_access = entry.access not in (Access.WRITE_ACCESS, Access.OWNER_ACCESS)

        if action == self.REVERT and no_access and not self.user.profile.is_moderator:
            raise forms.ValidationError(msg)

        if action == self.APPROVE and not self.user.profile.is_moderator:
            msg = "You have to be a moderator."
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


class PasteForm(forms.Form):
    """
    Used to paste analysis and files from clipboards
    """

    PASTE, CLEAR = 'PASTE', 'CLEAR'

    action = forms.CharField()

    def __init__(self, project, request, board, *args, **kwargs):

        self.project = project
        self.request = request
        self.board = board
        super().__init__(*args, **kwargs)


    def save(self):
        instance = self.cleaned_data.get("instance")
        action = self.cleaned_data.get('action')

        pasting_file = (isinstance(instance, Data) or isinstance(instance, Job))

        if action == self.CLEAR:
            self.request.session[self.board] = None

        if self.board == "recipe_clipboard" and action == self.PASTE:
            # Paste a recipe.
            attrs = auth.get_analysis_attr(instance, project=self.project)
            attrs.update(stream=instance.image, name=f"Copy of {instance.name}",
                         security=instance.security, user=self.request.user)
            new_recipe = auth.create_analysis(**attrs)
            # Ensure the diff gets inherited.
            new_recipe.last_valid = instance.last_valid
            new_recipe.save()
            messages.success(self.request, mark_safe(f"Pasted <b> {new_recipe.name}</b>"))

        if self.board == "files_clipboard" and pasting_file and action == self.PASTE:
            # Some files in clipboard might be outside job pathe
            files = self.request.session.get(self.board)
            files = [f for f in files if f.startswith(instance.get_data_dir())]
            type = '' if not isinstance(instance, Data) else instance.type
            # Add data to project
            data = auth.create_data(project=self.project, files=files, user=self.request.user,
                                    name=f'Copy of {instance.name}', type=type)
            messages.success(self.request, mark_safe(f"Pasted <b> {data.name}</b>"))

    def clean(self):

        cleaned_data = super(PasteForm, self).clean()
        current = self.request.session.get(self.board)
        instance = None
        access = {self.PASTE: Access.WRITE_ACCESS, self.CLEAR:Access.READ_ACCESS}

        if self.board == 'files_clipboard':
            if not auth.validate_files_clipboard(request=self.request) :
                raise forms.ValidationError("Invalid object in clipboard, copy again.")
            else:
                # Set instance to later access in save().
                instance = Job.objects.filter(uid=current[-1]).first()
                instance = instance or Data.objects.filter(uid=current[-1]).first()

        if self.board == 'recipe_clipboard':
            instance = Analysis.objects.filter(uid=current).first()

        # Make sure user has access to project
        has_access = auth.check_obj_access(instance=self.project, user=self.request.user, request=self.request,
                                           access=access[cleaned_data['action']])
        if not has_access:
            raise forms.ValidationError("Paste to another project.")

        # Update cleaned_data with instance.
        self.cleaned_data.update(dict(instance=instance))



class FileCopyForm(forms.Form):
    """
    Used to save paths found in jobs/data into files_clipboard"
    """

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
        paths = [ p for p in paths if os.path.exists(join(self.root_dir, p)) ]

        # Override cleaned_data to later access in save()
        self.cleaned_data = list(map(join, [self.root_dir] * len(paths), paths))


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
