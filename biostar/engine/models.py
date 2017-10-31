import logging, hjson

import mistune
from django.contrib.auth.models import Group
from django.contrib.auth.models import User
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.urls import reverse
from django.utils import timezone

from . import util, settings
from .const import *

logger = logging.getLogger("engine")

# The maximum length in characters for a typical name and text field.
MAX_NAME_LEN = 256
MAX_TEXT_LEN = 10000
MAX_LOG_LEN = 20 * MAX_TEXT_LEN

def join(*args):
    return os.path.abspath(os.path.join(*args))

class Bunch(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def make_html(text):
    return mistune.markdown(text)


def get_datatype(file):
    return Data.FILE


def filter_by_type():
    return


def data_upload_path(instance, filename):
    # Name the data by the filename.
    pieces = os.path.basename(filename).split(".")
    # File may have multiple extensions
    exts = ".".join(pieces[1:]) or "data"
    dataname = f"data-{instance.uid}.{exts}"
    return join(instance.project.get_path(), f"{instance.data_dir}", dataname)

def image_path(instance, filename):
    # Name the data by the filename.
    name, ext = os.path.splitext(filename)
    # File may have multiple extensions
    imgname = f"image-{instance.uid}{ext}"
    return  f"images/{imgname}"


class ProjectObjectManager(models.Manager):

    def get_queryset(self, user=Bunch(is_superuser=False)):

        if user.is_superuser:
            return super(ProjectObjectManager, self).get_queryset()

        return super(ProjectObjectManager, self).get_queryset().exclude(type=Project.ADMIN)


class ProjectAdminManager(models.Manager):

    def get_queryset(self):

        return super(ProjectAdminManager, self).get_queryset().filter(type=Project.ADMIN)


class Project(models.Model):
    ADMIN, USER = 1, 2
    PUBLIC, SHAREABLE, PRIVATE = 1,2,3
    PRIVACY_CHOICES = [ (PRIVATE, "Private"), (SHAREABLE, "Shareable Link"),  (PUBLIC, "Public") ]
    TYPE_CHOICES = [(ADMIN, "admin"), (USER, "user")]

    type = models.IntegerField(default=USER, choices=TYPE_CHOICES)
    privacy = models.IntegerField(default=SHAREABLE, choices=PRIVACY_CHOICES)

    image = models.ImageField(default=None, blank=True, upload_to=image_path)
    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')

    owner = models.ForeignKey(User)
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)

    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    # Project restircted to one group
    group = models.ForeignKey(Group)
    uid = models.CharField(max_length=32)

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

    # Override managers.
    objects = ProjectObjectManager()
    admins = ProjectAdminManager()

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.date = self.date or now
        self.html = make_html(self.text)

        # Takes first user group for now
        self.group = self.owner.groups.first()

        self.uid = self.uid or util.get_uuid(8)
        if not os.path.isdir(self.get_path()):
            os.makedirs(self.get_path())

        super(Project, self).save(*args, **kwargs)

    def __str__(self):
        return self.name

    def url(self):
        return reverse("project_view", kwargs=dict(id=self.id))

    def get_path(self):
        return join(settings.MEDIA_ROOT, "projects", f"proj-{self.uid}")


class Data(models.Model):
    ADMIN, USER = 1, 2
    FILE, COLLECTION = 1, 2
    PENDING, READY = 1, 2

    FILETYPE_CHOICES = [(FILE, "File"), (COLLECTION, "Collection")]
    TYPE_CHOICES = [(ADMIN, "Admin"), (USER, "User")]

    STATE_CHOICES = [(PENDING, "Pending"), (READY, "Ready")]

    type = models.IntegerField(default=USER, choices=TYPE_CHOICES)
    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')
    image = models.ImageField(default=None, blank=True, upload_to=image_path)

    owner = models.ForeignKey(User)
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)
    file_type = models.IntegerField(default=FILE, choices=FILETYPE_CHOICES)

    data_type = models.IntegerField(default=GENERIC_TYPE)
    project = models.ForeignKey(Project)
    size = models.CharField(null=True, max_length=256)

    state = models.IntegerField(default=PENDING, choices=STATE_CHOICES)
    file = models.FileField(null=True, upload_to=data_upload_path, max_length=500)
    uid = models.CharField(max_length=32)

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

    # Data directory.
    data_dir = models.FilePathField(default="")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.uid = self.uid or util.get_uuid(8)
        self.date = self.date or now
        self.html = make_html(self.text)

        # Build the data directory.
        data_dir = self.get_datadir()
        if not os.path.isdir(data_dir):
            os.makedirs(data_dir)
        self.data_dir = data_dir
        super(Data, self).save(*args, **kwargs)

    def peek(self):
        """
        Returns a preview of the data
        """
        return util.smart_preview(self.get_path())

    def set_size(self):
        """
        Sets the size of the data.
        """
        try:
            size = os.path.getsize(self.get_path())
        except:
            size = 0
        Data.objects.filter(id=self.id).update(size=size)

    def set_ready(self):
        Data.objects.filter(id=self.id).update(state=self.READY)

    def __str__(self):
        return self.name

    def get_datadir(self):
        return join(self.project.get_path(), f"store-{self.uid}")

    def get_path(self):
        return self.file.path

    def fill_dict(self, obj):
        """
        Mutates a dictionary to add more information.
        """

        # This is where the path is filled.
        obj['path'] = self.get_path()
        obj['data_id'] = self.id
        obj['name'] = self.name
        obj['uid'] = self.uid

class Analysis(models.Model):
    ADMIN, USER = 1, 2
    AUTHORIZED, UNDER_REVIEW = 1,2

    TYPE_CHOICES = [(ADMIN, "admin"), (USER, "user")]
    AUTH = [(AUTHORIZED,"authorized"), (UNDER_REVIEW, "under review")]

    type = models.IntegerField(default=USER, choices=TYPE_CHOICES)

    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    owner = models.ForeignKey(User)

    auth = models.IntegerField(default=UNDER_REVIEW, choices=AUTH)
    project = models.ForeignKey(Project)

    json_text = models.TextField(default="{}")
    template = models.TextField(default="makefile")
    uid = models.CharField(max_length=32)

    date = models.DateTimeField(auto_now_add=True, blank=True)
    image = models.ImageField(default=None, blank=True, upload_to=image_path)
    # Job start and end.
    start_date = models.DateTimeField(null=True, blank=True)
    end_date = models.DateTimeField(null=True, blank=True)

    # Will be false if the object is deleted.
    valid = models.BooleanField(default=True)

    def __str__(self):
        return self.name

    @property
    def json_data(self):
        "Returns the json_text as parsed json_data"
        return hjson.loads(self.json_text)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.uid = self.uid or util.get_uuid(8)
        self.date = self.date or now

        self.html = make_html(self.text)
        super(Analysis, self).save(*args, **kwargs)


class Job(models.Model):
    ADMIN, USER = 1, 2
    AUTHORIZED, UNDER_REVIEW = 1,2
    QUEUED, RUNNING, FINISHED, ERROR = 1, 2, 3, 4
    TYPE_CHOICES = [(ADMIN, "admin"), (USER, "user")]
    AUTH = [(AUTHORIZED, "authorized"), (UNDER_REVIEW, "under review")]
    STATE_CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"),
                     (FINISHED, "Finished"), (ERROR, "Error")]

    type = models.IntegerField(default=USER, choices=TYPE_CHOICES)
    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')
    image = models.ImageField(default=None, blank=True, upload_to=image_path)

    owner = models.ForeignKey(User)
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    analysis = models.ForeignKey(Analysis)
    project = models.ForeignKey(Project)
    json_text = models.TextField(default="commands")

    uid = models.CharField(max_length=32)
    template = models.TextField(default="makefile")
    auth = models.IntegerField(default=UNDER_REVIEW, choices=AUTH)

    # This will be set when the job attempts to run.
    script = models.TextField(default="")

    # Keeps track of errors.
    stdout_log = models.TextField(default="", max_length=MAX_LOG_LEN)

    # Standard error.
    stderr_log = models.TextField(default="", max_length=MAX_LOG_LEN)

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

    state = models.IntegerField(default=1, choices=STATE_CHOICES)

    path = models.FilePathField(default="")

    def is_running(self):
        return self.state == Job.RUNNING

    def __str__(self):
        return self.name

    def get_url(self, path=''):
        "Return the url to the job directory"
        return f"jobs/job-{self.uid}/" + path

    @property
    def json_data(self):
        "Returns the json_text as parsed json_data"
        return hjson.loads(self.json_text)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.name = self.name or self.analysis.name
        self.date = self.date or now
        self.html = make_html(self.text)

        self.uid = self.uid or util.get_uuid(8)
        self.template = self.analysis.template
        self.auth = self.analysis.auth

        self.name = self.name or self.analysis.name
        # write an index.html to the file
        if not os.path.isdir(self.path):
            path = join(settings.MEDIA_ROOT, "jobs", f"job-{self.uid}")
            os.makedirs(path)
            self.path = path

        super(Job, self).save(*args, **kwargs)

    def url(self):
        return reverse("job_view", kwargs=dict(id=self.id))


class Profile(models.Model):
    user = models.ForeignKey(User)


@receiver(post_save, sender=User)
def create_profile(sender, instance, created, **kwargs):
    if created:
        # Create a profile for user
        Profile.objects.create(user=instance)

        # Add every user to "public group"
        # instance.groups.add(Group.objects.get(name='Public'))


post_save.connect(create_profile, sender=User)
