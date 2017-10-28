import hjson
import logging

import mistune
from django.contrib.auth.models import Group
from django.contrib.auth.models import User
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.urls import reverse
from django.utils import timezone

from . import auth
from engine import settings
from . import util
from .const import *

logger = logging.getLogger("engine")

# The maximum length in characters for a typical name and text field.
MAX_NAME_LEN = 256
MAX_TEXT_LEN = 10000


def join(*args):
    return os.path.abspath(os.path.join(*args))


def make_html(text):
    return mistune.markdown(text)


def get_datatype(file):
    return Data.FILE


def filter_by_type():
    return


def upload_path(instance, filename):
    # Name the data by the filename.
    pieces = os.path.basename(filename).split(".")
    # File may have multiple extensions
    exts = ".".join(pieces[1:]) or "data"
    dataname = f"data-{instance.uid}.{exts}"
    return join(instance.project.get_path(), f"{instance.data_dir}", dataname)



class ProjectAdminManager(models.Manager):
    def get_queryset(self):
        return super(ProjectAdminManager, self).get_queryset().filter(type=Project.ADMIN)


class ProjectObjectManager(models.Manager):
    def get_queryset(self):
        return super(ProjectObjectManager, self).get_queryset().exclude(type=Project.ADMIN)


class Project(models.Model):
    ADMIN, USER = 1, 2
    TYPE_CHOICES = [(ADMIN, "admin"), (USER, "user")]
    type = models.IntegerField(default=USER, choices=TYPE_CHOICES)

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

    def get_data(self, data_type=None):
        """
        Returns a dictionary keyed by data stored in the project.
        """

        return auth.get_data(user=self.owner, project=self, data_type=data_type,
                             query=Data.objects.filter(project=self))

    def create_analysis(self, json_text, template, owner=None, summary='', name='', text='', type=None):
        """
        Creates analysis from a spec and template
        """
        analysis = auth.create_analysis(owner, self, Analysis, json_text, template, summary, name, text, type)
        logger.info(f"Created analysis id={analysis.id} of type_type={dict(Analysis.TYPE_CHOICES)[analysis.type]}")
        return analysis

    def create_data(self, stream=None, fname=None, name="data.bin", owner=None, text='', data_type=None, type=None):
        """
        Creates a data for the project from filename or a stream.
        """

        data = auth.create_data(owner, Data, self, stream, fname, name, text, data_type, type)

        logger.info(f"Added data id={data.id} of type={data.data_type}")
        return data


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

    owner = models.ForeignKey(User)
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)
    file_type = models.IntegerField(default=FILE, choices=FILETYPE_CHOICES)

    data_type = models.IntegerField(default=GENERIC_TYPE)
    project = models.ForeignKey(Project)
    size = models.CharField(null=True, max_length=256)

    state = models.IntegerField(default=PENDING, choices=STATE_CHOICES)
    file = models.FileField(null=True, upload_to=upload_path, max_length=500)
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
        obj['value'] = self.id
        obj['name'] = self.name


class Analysis(models.Model):
    ADMIN, USER = 1, 2
    TYPE_CHOICES = [(ADMIN, "admin"), (USER, "user")]
    type = models.IntegerField(default=USER, choices=TYPE_CHOICES)

    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    owner = models.ForeignKey(User)

    project = models.ForeignKey(Project)

    json_text = models.TextField(default="{}")
    template = models.TextField(default="makefile")
    uid = models.CharField(max_length=32)

    date = models.DateTimeField(auto_now_add=True, blank=True)

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

    def create_job(self, json_text='', json_data={}, owner=None, name=None, state=None, type=None):
        """
        Creates a job from an analysis.
        """
        job = auth.create_job(owner, self, Job, self.project, json_text, json_data, name, state, type)
        logger.info(f"Queued job: '{job.name}'")
        return job


class Job(models.Model):
    ADMIN, USER = 1, 2
    QUEUED, RUNNING, FINISHED, ERROR = 1, 2, 3, 4
    TYPE_CHOICES = [(ADMIN, "admin"), (USER, "user")]
    STATE_CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"),
                     (FINISHED, "Finished"), (ERROR, "Error")]

    type = models.IntegerField(default=USER, choices=TYPE_CHOICES)
    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')

    owner = models.ForeignKey(User)
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    analysis = models.ForeignKey(Analysis)
    project = models.ForeignKey(Project)
    json_text = models.TextField(default="commands")

    uid = models.CharField(max_length=32)
    template = models.TextField(default="makefile")

    # This will be set when the job attempts to run.
    script = models.TextField(default="")

    # Keeps track of errors.
    log = models.TextField(default="No data logged for current job", max_length=10 * MAX_TEXT_LEN)

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
        self.date = self.date or now
        self.html = make_html(self.text)

        self.uid = self.uid or util.get_uuid(8)
        self.template = self.analysis.template

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
