import copy
import shlex
import hjson
import io
import re

from django import forms
from django.template import Template, Context
from django.db.models import Sum
from django.utils.safestring import mark_safe
from django.utils.timezone import now
from django.contrib import messages
from django.urls import reverse
from django.core.validators import validate_slug

from django.conf import settings
from snowpenguin.django.recaptcha2.fields import ReCaptchaField
from snowpenguin.django.recaptcha2.widgets import ReCaptchaWidget
from biostar.accounts.models import User, Profile
from . import models, auth, factory, util
from .const import *
from .models import Project, Data, Analysis, Job, Access

# Share the logger with models.
logger = models.logger

TEXT_UPLOAD_MAX = 10000

# Maximum file size that can be uploaded to recipe in megabytes.
MAX_RECIPE_FILE_MB = 5


def join(*args):
    return os.path.abspath(os.path.join(*args))


def check_size(fobj, maxsize=0.3, field=None):
    # maxsize in megabytes!
    error_msg = ''
    try:
        if fobj and fobj.size > maxsize * 1024 * 1024.0:
            curr_size = fobj.size / 1024 / 1024.0
            prefix = f'{field} field : '.capitalize() if field else ''
            error_msg = prefix + f"file too large, {curr_size:0.1f}MB should be < {maxsize:0.1f}MB"

    except Exception as exc:
        error_msg = f"File size validation error: {exc}"

    if error_msg:
        raise forms.ValidationError(error_msg)

    return fobj


def check_upload_limit(file, user):
    """
    Checks if the file pushes user over their upload limit."
    """

    # Existing data.
    data = Data.objects.filter(owner=user, method=Data.UPLOAD)

    # The current cumulative size of the current data.
    current_size = data.aggregate(Sum("size"))["size__sum"] or 0

    # The projected size in MB.
    projected_size = file.size + current_size

    # Maximal cumulative sizes.
    max_size = user.profile.max_upload_size * 1024 * 1024

    # Allow much higher limits for staff.
    if user.is_staff:
        max_size = max_size * 10

    # Current file size in MB
    file_mb = file.size / 1024 / 1024

    # Verify projected data sizes.
    if projected_size > max_size:
        msg = f"You don't have enough storage space for data of size <b>{file_mb:.2f} MB</b>"
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


def add_captcha_field(request, fields):
    """Used to dynamically load captcha field into forms"""

    # Trusted users do not need a captcha check
    if request.user.is_authenticated and request.user.profile.trusted:
        return
    # Mutates the fields dict to add captcha field.
    if settings.RECAPTCHA_PRIVATE_KEY:
        fields["captcha"] = ReCaptchaField(widget=ReCaptchaWidget())
    return


class ProjectForm(forms.ModelForm):
    image = forms.ImageField(required=False)

    # Should not edit uid because data directories get recreated

    def __init__(self, request, create=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.request = request
        self.create = create

        # choices = filter(lambda pri: pri[0] != Project.SHAREABLE, Project.PRIVACY_CHOICES)

        self.fields["privacy"] = forms.CharField(widget=forms.Select(choices=Project.PRIVACY_CHOICES),
                                                 initial=self.instance.privacy,
                                                 required=False)

    class Meta:
        model = Project
        fields = ['name', 'text', 'privacy', 'rank', 'image']

    def clean_image(self):
        cleaned_data = super(ProjectForm, self).clean()
        image = cleaned_data.get('image')
        check_size(fobj=image)

        return image

    def clean(self):
        cleaned_data = super(ProjectForm, self).clean()

        user = self.request.user
        projects = Project.objects.filter(owner=user)
        privacy = cleaned_data.get("privacy") or 0

        # Trusted users are allowed everything.
        if user.is_authenticated and (user.is_staff or user.profile.trusted):
            return cleaned_data

        # Check project limit.
        if self.create and projects.count() > settings.MAX_PROJECTS:
            raise forms.ValidationError(
                f"You have exceeded the maximum number of projects allowed:{settings.MAX_PROJECTS}.")

        # Check privacy.
        if int(privacy) == Project.PUBLIC:
            raise forms.ValidationError(f"Only staff members can make public projects for now.")

        return cleaned_data

    def custom_save(self, owner):
        """Used to save on creation using custom function."""

        name = self.cleaned_data["name"]
        text = self.cleaned_data["text"]
        stream = self.cleaned_data["image"]
        project = auth.create_project(user=owner, name=name, text=text, stream=stream)
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
        type = self.cleaned_data["type"]
        name = self.cleaned_data['data_name']

        if stream:
            name = name or stream.name
        else:
            stream = io.StringIO(initial_value=input_text)

        data = auth.create_data(stream=stream, name=name, text=text, user=self.user,
                                project=self.project, type=type)
        if input_text and not self.cleaned_data["file"]:
            Data.objects.filter(pk=data.pk).update(method=Data.TEXTAREA)
            stream.close()

        return data

    class Meta:
        model = Data
        fields = ['data_name', 'input_text', 'text', "type"]

    def clean(self):

        cleaned_data = super(DataUploadForm, self).clean()
        upload = cleaned_data.get("file")
        text = cleaned_data.get("input_text")

        if not (upload or text):
            raise forms.ValidationError("Upload a file or write into the text field to create some data.")

        if upload:
            clean_file(fobj=upload, user=self.user,
                       project=self.project, check_name=False)

        else:
            if not cleaned_data.get("data_name"):
                raise forms.ValidationError("Name is required with text inputs.")

        total_count = Data.objects.filter(owner=self.user).count()
        if total_count >= settings.MAX_DATA:
            raise forms.ValidationError(f"Exceeded maximum amount of data:{settings.MAX_DATA}.")
        return cleaned_data

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

        self.instance.lastedit_user = self.user
        self.instance.lasedit_date = now()
        Project.objects.filter(uid=self.instance.project.uid).update(lastedit_user=self.user,
                                                                     lastedit_date=now()
                                                                     )

        return super(DataEditForm, self).save(commit)

    class Meta:
        model = Data
        fields = ['name', 'text', "type"]

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


class RecipeForm(forms.ModelForm):
    """
    Fields that are not submitted are set to existing values.
    """
    image = forms.ImageField(required=False)
    uid = forms.CharField(max_length=32, validators=[validate_slug], required=True)
    json_text = forms.CharField(max_length=MAX_TEXT_LEN, initial="{}", required=True)
    template = forms.CharField(max_length=MAX_TEXT_LEN, initial="# Code goes here", required=True)
    name = forms.CharField(max_length=MAX_NAME_LEN, required=True)
    rank = forms.FloatField(required=True, initial=100)
    text = forms.CharField(initial="Recipe description", widget=forms.Textarea, required=True)

    def __init__(self, user, creating=False, project=None, *args, **kwargs):
        self.creating = creating
        self.user = user
        self.project = project

        super().__init__(*args, **kwargs)

        # Admins get an added field
        if self.user.is_superuser:
            authorized = self.instance.runnable()
            self.fields['authorized'] = forms.BooleanField(widget=forms.CheckboxInput(attrs={'class': 'ui checkbox'}),
                                                           initial=authorized, required=False)

    class Meta:
        model = Analysis
        fields = ["name", "rank", "text", "uid", "json_text", "template"]

    def get_initial(self):
        """
        Returns the initial data to use for forms on this view.
        """
        initial = super(RecipeForm, self).get_initial()
        for field in self.Meta.fields:
            initial['field'] = getattr(self.instance, field)

        return initial

    def validate_readable(self):

        is_readable = auth.is_readable(user=self.user, project=self.project)
        if not is_readable:
            raise forms.ValidationError('You need read access to the project create a recipe.')

    def validate_writable(self):
        # Check write access when editing
        is_writable = auth.writeable_recipe(user=self.user, source=self.instance, project=self.project)
        if not is_writable:
            raise forms.ValidationError('You need write access to the original recipe to edit.')

    def clean(self):
        """
        Applies security measures to recipe editing.
        """
        cleaned_data = super(RecipeForm, self).clean()

        if self.user.is_anonymous:
            raise forms.ValidationError('You need to be logged in.')

        if self.creating:
            # Check if recipe is readable to user
            self.validate_readable()
        else:
            # Check to see if the
            self.validate_writable()

        # Fill in not submitted fields.
        for field in self.Meta.fields:
            cleaned_data['field'] = getattr(self.instance, field)

        template = self.cleaned_data.get('template') or self.instance.template
        json_text = self.cleaned_data.get('json_text') or self.instance.json_text

        # Shortcuts to security conditions.
        template_changed = (template != self.instance.template)
        json_changed = (json_text != self.instance.json_text)

        not_superuser = not self.user.is_superuser

        # Update the recipe security when template or JSON have been
        # touched by non admin users.
        if (json_changed or template_changed) and not_superuser and not self.creating:
            Analysis.objects.filter(uid=self.instance.uid).update(security=Analysis.NOT_AUTHORIZED)
            self.instance.security = Analysis.NOT_AUTHORIZED
            self.instance.save()

        return cleaned_data

    def clean_image(self):
        cleaned_data = super(RecipeForm, self).clean()
        image = cleaned_data.get('image')
        check_size(fobj=image)

        return image

    def clean_json_text(self):
        cleaned_data = super(RecipeForm, self).clean()
        json_text = cleaned_data.get('json_text')

        # Ensure correct JSON syntax.
        try:
            hjson.loads(json_text)
        except Exception as exc:
            raise forms.ValidationError(f'Error with recipe JSON:{exc}')

        return json_text

    def clean_uid(self):
        cleaned_data = super(RecipeForm, self).clean()
        uid = cleaned_data.get('uid')

        if self.creating:
            # Check if uid already exists when creating a recipe.
            recipe = Analysis.objects.filter(uid=uid).first()
            if recipe:
                raise forms.ValidationError("Recipe uid already exists.")
        return uid

    def save(self, commit=True):

        self.instance.lastedit_date = now()
        self.instance.lastedit_user = self.user
        Project.objects.filter(uid=self.instance.project.uid).update(lastedit_date=now(),
                                                                     lastedit_user=self.user)
        image = self.cleaned_data['image']
        self.instance.image = image or self.instance.image
        if self.user.is_superuser:
            security = Analysis.AUTHORIZED if self.cleaned_data.get('authorized') else Analysis.NOT_AUTHORIZED
            self.instance.security = security

        return super().save(commit)


class JobEditForm(forms.ModelForm):

    def __init__(self, user, *args, **kwargs):
        self.user = user
        super().__init__(*args, **kwargs)

    class Meta:
        model = Job
        fields = ['name', "image", 'text']

    def save(self, commit=True):
        self.instance.lastedit_date = now()
        self.instance.lastedit_user = self.user
        Project.objects.filter(uid=self.instance.project.uid).update(lastedit_date=now(),
                                                                     lastedit_user=self.user)
        return super().save(commit)


def clean_text(textbox):
    return shlex.quote(textbox)


class RecipeInterface(forms.Form):
    # The name of results when running the recipe.
    # name = forms.CharField(max_length=256, label="Name", help_text="This is how you can identify the run.")

    def __init__(self, request, json_data, analysis=None, project=None, add_captcha=True, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # The json data determines what fields does the form have.
        self.json_data = json_data

        # The project is required to select data from.
        self.analysis = analysis
        self.project = analysis.project if analysis else project

        # Get request specific information
        self.request = request
        self.user = self.request.user

        # Create the dynamic field from each key in the data.
        for name, data in self.json_data.items():
            field = factory.dynamic_field(data, self.project)

            # Insert only valid fields.
            if field:
                self.fields[name] = field

        if 0:
            add_captcha_field(request=request, fields=self.fields)

    def clean(self):

        # Validate default fields.
        super(RecipeInterface, self).clean()
        valid, msg = auth.validate_recipe_run(user=self.user, recipe=self.analysis)
        if not valid:
            raise forms.ValidationError(msg)

        for field, item in self.json_data.items():

            if field in self.cleaned_data and item.get('display') == UPLOAD:
                stream = self.request.FILES.get(field)
                if not stream:
                    continue
                # Validate the file size.
                check_size(stream, field=field, maxsize=MAX_RECIPE_FILE_MB)

        self.validate_text_fields()

    def validate_text_fields(self):
        """Validate Character fields """
        bool_map = {'false': False, 'true': True}
        # Default pattern matches any alphanumeric string with a given length
        default_pattern = r"^\w{1,20}$"

        for field, item in self.json_data.items():
            val = self.cleaned_data.get(field)

            # Validate text fields
            if (val is None) or (item.get("display") != TEXTBOX):
                continue

            # Acceptable regex pattern
            regex_pattern = item.get("regex", default_pattern)

            if re.fullmatch(regex_pattern, val) is None:
                msg = f"{field} : contains invalid patterns. Valid pattern:{regex_pattern}."
                raise forms.ValidationError(msg)