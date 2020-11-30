import logging

import toml as hjson
import mistune
import urllib.parse
import base64
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.template import loader
from django.urls import reverse
from django.utils import timezone
from django.conf import settings
from biostar.accounts.models import User
from . import util
from .const import *

logger = logging.getLogger("engine")


def join(*args):
    return os.path.abspath(os.path.join(*args))


class Bunch(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def make_html(text, user=None):
    if user and user.profile.trusted:
        html = mistune.markdown(text, escape=False)
    else:
        html = mistune.markdown(text, escape=True)
    return html


def image_path(instance, filename):
    # Name the data by the filename.
    name, ext = os.path.splitext(filename)

    uid = util.get_uuid(6)
    dirpath = instance.get_project_dir()
    imgname = f"images/image-{uid}{ext}"

    # Uploads need to go relative to media directory.
    path = os.path.relpath(dirpath, settings.MEDIA_ROOT)

    imgpath = os.path.join(path, imgname)

    return imgpath


def snippet_images(instance, filename):
    # Name the data by the filename.
    name, ext = os.path.splitext(filename)
    uid = instance.uid or util.get_uuid(6)
    imgname = f"image-{uid}{ext}"

    dirpath = os.path.abspath(os.path.join(settings.MEDIA_ROOT, 'snippets', 'images'))

    # Uploads need to go relative to media directory.
    path = os.path.relpath(dirpath, settings.MEDIA_ROOT)
    imgpath = os.path.join(path, imgname)

    return imgpath


class Manager(models.Manager):

    def get_queryset(self):
        "Regular queries exclude deleted stuff"
        return super().get_queryset().select_related("owner", "owner__profile", "lastedit_user",
                                                     "lastedit_user__profile")


class SnippetType(models.Model):
    image = models.ImageField(default=None, blank=True, upload_to=snippet_images, max_length=MAX_FIELD_LEN)
    uid = models.CharField(max_length=MAX_TEXT_LEN, unique=True)
    # The name referring to this
    name = models.CharField(max_length=MAX_NAME_LEN)
    owner = models.ForeignKey(User, on_delete=models.CASCADE)

    # Appears to all uses
    default = models.BooleanField(default=False)

    def save(self, *args, **kwargs):
        self.uid = self.uid or util.get_uuid(6)
        super(SnippetType, self).save(*args, **kwargs)


class Snippet(models.Model):
    help_text = models.CharField(max_length=MAX_TEXT_LEN)
    uid = models.CharField(max_length=MAX_TEXT_LEN, unique=True)
    owner = models.ForeignKey(User, on_delete=models.CASCADE)
    command = models.CharField(max_length=MAX_TEXT_LEN, null=True)
    # Link a command to one type
    type = models.ForeignKey(SnippetType, on_delete=models.CASCADE)

    # Appears to all uses
    default = models.BooleanField(default=False)

    def save(self, *args, **kwargs):
        self.uid = self.uid or util.get_uuid(6)
        self.owner = self.owner or self.type.owner
        super(Snippet, self).save(*args, **kwargs)


def img_to_str(img, decode=True):
    exists = img and os.path.exists(img.path)

    img = img.path if exists else os.path.join(settings.STATIC_ROOT, "images", "placeholder.png")
    stream = open(img, 'rb').read()
    # Convert image to base64 string
    strimg = base64.b64encode(stream)
    if decode:
        # Decode from bytes to string, enables JSON serialization.
        strimg = strimg.decode()
    return strimg


class Project(models.Model):
    PUBLIC, SHAREABLE, PRIVATE = 1, 2, 3
    PRIVACY_CHOICES = [(PRIVATE, "Private"), (SHAREABLE, "Shared"), (PUBLIC, "Public")]

    # Rank in a project list.
    rank = models.FloatField(default=100)

    # The user that edited the object most recently.
    lastedit_user = models.ForeignKey(User, related_name='proj_editor', null=True, on_delete=models.CASCADE)
    lastedit_date = models.DateTimeField(default=timezone.now)

    # Limits who can access the project.
    privacy = models.IntegerField(default=PRIVATE, choices=PRIVACY_CHOICES)
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN)
    name = models.CharField(default="New Project", max_length=MAX_NAME_LEN)
    deleted = models.BooleanField(default=False)

    # We need to keep the owner.
    owner = models.ForeignKey(User, null=False, on_delete=models.CASCADE)
    text = models.TextField(default='Project description.', max_length=MAX_TEXT_LEN)

    html = models.TextField(default='html', max_length=MAX_LOG_LEN)
    date = models.DateTimeField(auto_now_add=True)
    # Internal uid that is not editable.
    uid = models.CharField(max_length=32, unique=True)

    # FilePathField points to an existing project directory.
    dir = models.FilePathField(max_length=MAX_FIELD_LEN, default='')

    sharable_token = models.CharField(max_length=32, null=True, unique=True)

    data_count = models.IntegerField(default=0, null=True, db_index=True)
    recipes_count = models.IntegerField(default=0, null=True, db_index=True)
    jobs_count = models.IntegerField(default=0, null=True, db_index=True)

    objects = Manager()

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.date = self.date or now
        self.sharable_token = self.sharable_token or util.get_uuid(30)
        self.html = make_html(self.text, user=self.lastedit_user)
        self.name = self.name[:MAX_NAME_LEN]
        self.uid = self.uid or util.get_uuid(8)
        self.lastedit_user = self.lastedit_user or self.owner
        self.lastedit_date = now

        super(Project, self).save(*args, **kwargs)

    def __str__(self):
        return self.name

    def uid_is_set(self):
        assert bool(self.uid.strip()), "Sanity check. UID should always be set."

    def url(self):
        self.uid_is_set()
        return reverse("project_view", kwargs=dict(uid=self.uid))

    def get_project_dir(self):
        return self.dir

    def get_data_dir(self):
        "Match consistency of data dir calls"
        return self.get_project_dir()

    def set_counts(self, save=True):
        """
        Set the data, recipe, and job count for this project
        """

        # Set the counts on current instance and also update in database.
        self.recipes_count = self.analysis_set.filter(project=self, deleted=False).count()
        self.data_count = self.data_set.filter(deleted=False).count()
        self.job_count = self.job_set.filter(deleted=False).count()

        Project.objects.filter(id=self.id).update(data_count=self.data_count,
                                                  recipes_count=self.recipes_count,
                                                  jobs_count=self.job_count)

    @property
    def is_public(self):
        return self.privacy == self.PUBLIC

    @property
    def is_private(self):
        return self.privacy == self.PRIVATE

    @property
    def project(self):
        return self

    @property
    def json_text(self):
        return hjson.dumps(self.api_data)

    @property
    def api_data(self):
        img = self.image
        strimg = img_to_str(img=img)

        json_data = dict(
            uid=self.uid,
            name=self.name,
            privacy=dict(self.PRIVACY_CHOICES)[self.privacy],
            text=self.text,
            url=settings.BASE_URL,
            project_uid=self.uid,
            id=self.pk,
            image=strimg,
            # Insert recipes API data in there as well.
            recipes=[
                r.api_data for r in self.analysis_set.all()
            ]
        )

        return json_data

    @property
    def summary(self):
        """
        Returns first line of text
        """
        lines = self.text.splitlines() or ['']
        first = lines[0]
        return first

    @property
    def delete_url(self):
        return reverse('project_delete', kwargs=dict(uid=self.uid))

    @property
    def is_shareable(self):
        return self.privacy == self.SHAREABLE

    def get_sharable_link(self):

        # Return a sharable link if the project is shareable
        if self.is_shareable:
            return reverse('project_share', kwargs=dict(token=self.sharable_token))

        return '/'

    def get_name(self):
        if self.deleted:
            return f'Deleted Project: {self.name}'
        return self.name


class Access(models.Model):
    """
    Allows access of users to Projects.
    """
    NO_ACCESS, READ_ACCESS, WRITE_ACCESS, SHARE_ACCESS = 1, 2, 3, 4
    ACCESS_CHOICES = [
        (NO_ACCESS, "No Access"),
        (READ_ACCESS, "Read Access"),
        (WRITE_ACCESS, "Write Access"),
        (SHARE_ACCESS, "Share Access"),
    ]

    ACCESS_MAP = dict(ACCESS_CHOICES)

    user = models.ForeignKey(User, on_delete=models.CASCADE)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)

    access = models.IntegerField(default=NO_ACCESS, choices=ACCESS_CHOICES, db_index=True)
    date = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.user} on {self.project.name}"

    def save(self, *args, **kwargs):
        super(Access, self).save(*args, **kwargs)


@receiver(post_save, sender=Project)
def update_access(sender, instance, created, raw, update_fields, **kwargs):
    # Give the owner WRITE ACCESS if they do not have it.
    entry = Access.objects.filter(user=instance.owner, project=instance, access=Access.WRITE_ACCESS)
    if entry.first() is None:
        entry = Access.objects.create(user=instance.owner, project=instance, access=Access.WRITE_ACCESS)


class Data(models.Model):
    PENDING, READY, ERROR, = 1, 2, 3
    STATE_CHOICES = [(PENDING, "Pending"), (READY, "Ready"), (ERROR, "Error")]
    state = models.IntegerField(default=PENDING, choices=STATE_CHOICES)

    LINK, UPLOAD, TEXTAREA = 1, 2, 3
    METHOD_CHOICE = [(LINK, "Linked Data"), (UPLOAD, "Uploaded Data"), (TEXTAREA, "Text Field")]
    method = models.IntegerField(default=LINK, choices=METHOD_CHOICE)

    name = models.CharField(max_length=MAX_NAME_LEN, default="My Data")
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN)

    deleted = models.BooleanField(default=False)

    # Rank on a data list.
    rank = models.FloatField(default=100)

    # The user that edited the object most recently.
    lastedit_user = models.ForeignKey(User, related_name='data_editor', null=True, on_delete=models.CASCADE)
    lastedit_date = models.DateTimeField(default=timezone.now)

    owner = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    text = models.TextField(default='Data description.', max_length=MAX_TEXT_LEN, blank=True)
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    type = models.CharField(max_length=MAX_NAME_LEN, default="DATA")
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    size = models.BigIntegerField(default=0)

    # FilePathField points to an existing file path
    file = models.FilePathField(max_length=MAX_FIELD_LEN, path='')

    # FilePathField points to the data directory
    dir = models.FilePathField(max_length=MAX_FIELD_LEN, default='')

    # FilePathField points to the toc file.
    toc = models.FilePathField(max_length=MAX_FIELD_LEN, default='')

    # Get the file count from the toc file.
    file_count = models.IntegerField(default=0)

    uid = models.CharField(max_length=32, unique=True)

    objects = Manager()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.name = self.name[:MAX_NAME_LEN]
        self.uid = self.uid or util.get_uuid(8)
        self.date = self.date or now
        self.html = make_html(self.text, user=self.lastedit_user)
        self.owner = self.owner or self.project.owner
        self.type = self.type.replace(" ", '')
        self.lastedit_user = self.lastedit_user or self.owner or self.project.owner
        self.lastedit_date = now

        super(Data, self).save(*args, **kwargs)

        # Set the counts
        self.project.set_counts(save=True)

    def peek(self):
        """
        Returns a preview of the data
        """
        try:
            target = self.get_path()
            lines = open(target, 'rt').readlines()
            if len(lines) == 1:
                target = lines[0]
                return util.smart_preview(target)
            else:
                data_dir = self.get_data_dir()
                rels = [os.path.relpath(path, data_dir) for path in lines]
                return "".join(rels)
        except Exception as exc:
            return f"Error :{exc}"

    def table_of_contents(self):
        try:
            target = self.get_path()
            lines = open(target, 'rt').readlines()
            lines = [os.path.relpath(path, self.get_data_dir()) for path in lines]
        except Exception as exc:
            return f"Error :{exc}"

        return lines

    def __str__(self):
        return self.name

    @property
    def delete_url(self):
        return reverse('data_delete', kwargs=dict(uid=self.uid))

    def get_data_dir(self):
        "The data directory"
        return self.dir

    def get_project_dir(self):
        return self.project.get_project_dir()

    def get_path(self):
        return self.toc

    def make_toc(self):

        tocname = self.get_path()

        collect = util.findfiles(self.get_data_dir(), collect=[])

        # Create a sorted file path collection.
        collect.sort()
        # Write the table of contents.
        with open(tocname, 'w') as fp:
            fp.write("\n".join(collect))

        # Find the cumulative size of the files.
        size = 0
        for elem in collect:
            if os.path.isfile(elem):
                size += os.stat(elem, follow_symlinks=True).st_size

        self.size = size
        self.file = tocname
        self.file_count = len(collect)
        Data.objects.filter(id=self.id).update(size=self.size, file=self.file, file_count=self.file_count)

        return tocname

    def can_unpack(self):
        cond = str(self.get_path()).endswith("tar.gz")
        return cond

    def get_files(self):
        fnames = [line.strip() for line in open(self.get_path(), 'rt')]
        return fnames if len(fnames) else [""]

    def get_url(self, path=""):
        "Returns url to the data directory"
        return f"projects/{self.project.uid}/{self.uid}/" + path

    def url(self):
        return reverse('data_view', kwargs=dict(uid=self.uid))

    def fill_dict(self, obj):
        """
        Mutates a dictionary object to fill in more fields based
        on the current object.
        """
        fnames = self.get_files()
        if fnames:
            obj['value'] = fnames[0]
        else:
            obj['value'] = 'MISSING'

        obj['files'] = fnames
        obj['toc'] = self.get_path()
        obj['file_list'] = self.get_path()
        obj['id'] = self.id
        obj['name'] = self.name
        obj['uid'] = self.uid
        obj['data_dir'] = self.get_data_dir()
        obj['project_dir'] = self.get_project_dir()
        obj['data_url'] = self.url()

    @property
    def summary(self):
        """
        Returns first line of text
        """
        lines = self.text.splitlines() or ['']
        first = lines[0]
        return first

    def get_name(self):
        if self.deleted:
            return f'Deleted Data: {self.name}'

        return self.name


class Analysis(models.Model):
    AUTHORIZED, NOT_AUTHORIZED = 1, 2

    SECURITY_STATES = [
        (AUTHORIZED, "Trusted users"),
        (NOT_AUTHORIZED, "Admin only")
    ]

    security = models.IntegerField(default=NOT_AUTHORIZED, choices=SECURITY_STATES)

    deleted = models.BooleanField(default=False)
    uid = models.CharField(max_length=32, unique=True)

    name = models.CharField(max_length=MAX_NAME_LEN, default="New Recipe")
    text = models.TextField(default='This is the recipe description.', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    owner = models.ForeignKey(User, on_delete=models.CASCADE)

    # The rank in a recipe list.
    rank = models.FloatField(default=100)

    # Root recipe this recipe has been copied from
    root = models.ForeignKey('self', null=True, blank=True, on_delete=models.SET_NULL)

    # The user that edited the object most recently.
    lastedit_user = models.ForeignKey(User, related_name='analysis_editor', null=True, on_delete=models.CASCADE)
    lastedit_date = models.DateTimeField(default=timezone.now)

    project = models.ForeignKey(Project, on_delete=models.CASCADE)

    json_text = models.TextField(default="", max_length=MAX_TEXT_LEN)

    template = models.TextField(default="")
    last_valid = models.TextField(default='')

    date = models.DateTimeField(auto_now_add=True, blank=True)
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN,
                              help_text="Optional image")

    objects = Manager()

    def __str__(self):
        return self.name

    @property
    def json_data(self):
        """
        Returns the json_text as parsed json_data
        """
        try:
            json_data = hjson.loads(self.json_text)
        except Exception as exc:
            logger.error(f"{exc}. json_text={self.json_text}")
            json_data = {}

        # Generates file names
        base = f"{'_'.join(self.name.split())}_{self.project.uid}_{self.pk}"
        template_name = f"{base}.sh"
        image_name = f"{base}.png"

        # Previously set values.
        current_settings = json_data.get("settings") or {}

        # Overwrite any previously set values with current information.
        current_settings["name"] = self.name
        current_settings["template"] = template_name
        current_settings["image"] = image_name
        current_settings["id"] = self.pk
        current_settings["recipe_uid"] = self.uid
        current_settings["uid"] = self.uid
        current_settings["help"] = self.text
        current_settings["url"] = settings.BASE_URL
        current_settings['root_id'] = self.root.id if self.root else ""
        current_settings['root_uid'] = self.root.uid if self.root else ""
        # Put them back into settings.
        json_data["settings"] = current_settings

        return json_data

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.uid = self.uid or util.get_uuid(8)
        self.date = self.date or now
        self.text = self.text or "Recipe description"
        self.name = self.name[:MAX_NAME_LEN] or "New Recipe"
        self.html = make_html(self.text, user=self.lastedit_user)
        self.lastedit_user = self.lastedit_user or self.owner or self.project.owner
        self.lastedit_date = now

        # Ensure Unix line endings.
        self.template = self.template.replace('\r\n', '\n') if self.template else ""

        self.project.set_counts(save=True)
        super(Analysis, self).save(*args, **kwargs)

    @property
    def api_data(self):
        img = self.image
        strimg = img_to_str(img=img)
        payload = {"uid": self.uid,
                   "name": self.name,
                   "text": self.text,
                   "json": hjson.dumps(self.json_data),
                   "template": self.template,
                   "image": strimg}

        return payload

    def get_project_dir(self):
        return self.project.get_project_dir()

    @property
    def delete_url(self):
        return reverse('recipe_delete', kwargs=dict(uid=self.uid))

    @property
    def is_cloned(self):
        """
        Return True if recipe is a clone ( linked ).
        """
        return self.root is not None

    @property
    def is_root(self):
        """
        Return True if recipe is a root.
        """
        return self.root is None

    @property
    def is_clone(self):
        """
        Return True if recipe is a clone.
        """
        return not self.is_root

    def update_children(self):
        """
        Update information for children belonging to this root.
        """
        # Get all children of this root
        children = Analysis.objects.filter(root=self)

        # Update all children information
        children.update(json_text=self.json_text,
                        template=self.template,
                        name=self.name,
                        security=self.security,
                        lastedit_date=self.lastedit_date,
                        lastedit_user=self.lastedit_user,
                        text=self.text,
                        html=self.html,
                        image=self.image)

        # Update last edit user and date for children projects.
        Project.objects.filter(analysis__root=self).update(lastedit_date=self.lastedit_date,
                                                           lastedit_user=self.lastedit_user)

    def url(self):
        assert self.uid, "Sanity check. UID should always be set."
        return reverse("recipe_view", kwargs=dict(uid=self.uid))

    def runnable(self):
        return self.security == self.AUTHORIZED

    def edit_url(self):
        # Return root edit url if this recipe is cloned.
        # if self.is_cloned:
        #    return reverse('recipe_edit', kwargs=dict(uid=self.root.uid))

        return reverse('recipe_view', kwargs=dict(uid=self.uid))

    @property
    def summary(self):
        """
        Returns first line of text
        """
        lines = self.text.splitlines() or ['']
        first = lines[0]
        return first

    def download_url(self):
        fname = "_".join(self.name.split()) + ".sh"
        return reverse('recipe_download', kwargs=dict(uid=self.uid, fname=fname))

    def get_name(self):
        if self.deleted:
            return f'Deleted Recipe: {self.name}'

        return self.name


class Job(models.Model):
    AUTHORIZED, UNDER_REVIEW = 1, 2
    AUTH_CHOICES = [(AUTHORIZED, "Authorized"), (UNDER_REVIEW, "Authorization Required")]

    QUEUED, RUNNING, COMPLETED, ERROR, SPOOLED, PAUSED = range(1, 7)

    STATE_CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"), (PAUSED, "Paused"),
                     (SPOOLED, "Spooled"), (COMPLETED, "Completed"), (ERROR, "Error")]

    state = models.IntegerField(default=QUEUED, choices=STATE_CHOICES)

    deleted = models.BooleanField(default=False)
    name = models.CharField(max_length=MAX_NAME_LEN, default="New results")
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN)

    # The user that edited the object most recently.
    lastedit_user = models.ForeignKey(User, related_name='job_editor', null=True, on_delete=models.CASCADE)
    lastedit_date = models.DateTimeField(default=timezone.now)

    owner = models.ForeignKey(User, on_delete=models.CASCADE)
    text = models.TextField(default='Results generated by running the recipe.', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')

    # Job creation date
    date = models.DateTimeField(auto_now_add=True)

    # Job runtime date.
    start_date = models.DateTimeField(null=True, blank=True)
    end_date = models.DateTimeField(null=True, blank=True)

    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    json_text = models.TextField(default="commands")

    uid = models.CharField(max_length=32)

    template = models.TextField(default="makefile")

    # Set the security level.
    security = models.IntegerField(default=UNDER_REVIEW, choices=AUTH_CHOICES)

    # This will be set when the job attempts to run.
    script = models.TextField(default="")

    # Keeps track of errors.
    stdout_log = models.TextField(default="", max_length=MAX_LOG_LEN)

    # Standard error.
    stderr_log = models.TextField(default="", max_length=MAX_LOG_LEN)

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

    path = models.FilePathField(default="")

    objects = Manager()

    def is_running(self):
        return self.state == Job.RUNNING

    def is_success(self):
        return self.state == Job.COMPLETED

    def is_error(self):
        return self.state == Job.ERROR

    def is_started(self):
        """
        This job has been initiated.
        """
        return self.state in [Job.QUEUED, Job.SPOOLED, Job.RUNNING]

    def is_finished(self):
        """
        This job is fishined
        """
        return self.state in [Job.ERROR, Job.COMPLETED]

    def __str__(self):
        return self.name

    @property
    def delete_url(self):
        return reverse('job_delete', kwargs=dict(uid=self.uid))

    def get_url(self, path=''):
        """
        Return the url to the job directory
        """
        return f"jobs/{self.uid}/" + path

    def url(self):
        url = reverse("job_view", kwargs=dict(uid=self.uid))

        # Anchor to the logs when the job is running
        if self.is_running():
            url = f"{url}#log/"

        return url

    def get_project_dir(self):
        return self.project.get_project_dir()

    def get_data_dir(self):
        return self.path

    @property
    def json_data(self):
        "Returns the json_text as parsed json_data"
        try:
            data_dict = hjson.loads(self.json_text)
        except Exception as exc:
            logger.error(f"{exc}; text={self.json_text}")
            data_dict = {}

        return data_dict

    def elapsed(self):
        if not (self.start_date and self.end_date):
            value = ''
        else:
            seconds = int((self.end_date - self.start_date).seconds)
            if seconds < 60:
                value = f'{seconds} seconds'
            elif seconds < 3600:
                minutes = int(seconds / 60)
                value = f'{minutes} minutes'
            else:
                hours = round(seconds / 3600, 1)
                value = f'{hours} hours'

        return value

    def done(self):
        return self.state == Job.COMPLETED

    def make_path(self):
        path = join(settings.MEDIA_ROOT, "jobs", f"{self.uid}")
        return path

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.name = self.name or f"Results for: {self.analysis.name}"
        self.date = self.date or now
        self.text = self.text or self.analysis.text
        self.html = make_html(self.text, user=self.lastedit_user)
        self.name = self.name[:MAX_NAME_LEN]
        self.uid = self.uid or util.get_uuid(8)
        self.template = self.analysis.template
        self.stderr_log = self.stderr_log[:MAX_LOG_LEN]
        self.stdout_log = self.stdout_log[:MAX_LOG_LEN]
        self.name = self.name or self.analysis.name

        self.lastedit_user = self.lastedit_user or self.owner or self.project.owner
        self.lastedit_date = now

        super(Job, self).save(*args, **kwargs)

    @property
    def parameter_summary(self):
        """
        """
        summary_template = "widgets/job_summary.md"
        context = dict(data=self.json_data)
        template = loader.get_template(summary_template)
        result = template.render(context)
        return result

    @property
    def summary(self):
        """
        Return the first line in self.text.
        """
        lines = self.text.splitlines() or ['']
        first = lines[0]
        return first

    def runnable(self):
        """
        Job is authorized to run
        """
        authorized = self.analysis.runnable() and self.security == self.AUTHORIZED
        return authorized

    def get_name(self):
        if self.deleted:
            return f'Deleted Job: {self.name}'

        return self.name
