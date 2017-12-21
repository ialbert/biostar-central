import logging

import hjson
import mistune
from django.db import models
from django.db.models.signals import post_save
from django.db.models.signals import pre_save
from django.dispatch import receiver
from django.urls import reverse
from django.utils import timezone

from biostar import settings
from biostar.tools import const
from biostar.accounts.models import User, Group
from . import util
from .const import *

logger = logging.getLogger("engine")

# The maximum length in characters for a typical name and text field.
MAX_NAME_LEN = 256
MAX_FIELD_LEN = 1024
MAX_TEXT_LEN = 10000
MAX_LOG_LEN = 20 * MAX_TEXT_LEN


def join(*args):
    return os.path.abspath(os.path.join(*args))


class Bunch(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def make_html(text):
    return mistune.markdown(text)


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


class Project(models.Model):
    PUBLIC, SHAREABLE, PRIVATE = 1, 2, 3
    PRIVACY_CHOICES = [(PRIVATE, "Private"), (SHAREABLE, "Shareable Link"), (PUBLIC, "Public")]

    ACTIVE, DELETED = 1, 2
    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted")]

    # Affects the sort order.
    sticky = models.BooleanField(default=False)

    # Limits who can access the project.
    privacy = models.IntegerField(default=PRIVATE, choices=PRIVACY_CHOICES)

    # The state a project may be in.
    state = models.IntegerField(default=ACTIVE, choices=STATE_CHOICES)

    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN)
    name = models.CharField(default="name", max_length=MAX_NAME_LEN)
    summary = models.TextField(default='summary', max_length=MAX_TEXT_LEN)

    # We need to keep the owner.
    owner = models.ForeignKey(User, null=False)
    text = models.TextField(default='description', max_length=MAX_TEXT_LEN)

    html = models.TextField(default='html', max_length=MAX_LOG_LEN)
    date = models.DateTimeField(auto_now_add=True)

    uid = models.CharField(max_length=32, unique=True)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.date = self.date or now
        self.html = make_html(self.text)
        self.name = self.name[:MAX_NAME_LEN]
        self.uid = self.uid or util.get_uuid(8)
        if not os.path.isdir(self.get_project_dir()):
            os.makedirs(self.get_project_dir())

        super(Project, self).save(*args, **kwargs)

    def __str__(self):
        return self.name

    def url(self):
        return reverse("project_view", kwargs=dict(uid=self.uid))

    def get_project_dir(self):
        return join(settings.MEDIA_ROOT, "projects", f"proj-{self.uid}")


class DataType(models.Model):

    name = models.CharField(default="name", max_length=MAX_NAME_LEN)

    # Symobol is what we enter in the json file
    symbol = models.CharField(max_length=MAX_FIELD_LEN)

    # given numeric value if one not given
    numeric = models.IntegerField(auto_created=True)

    help = models.CharField(default="description", max_length=MAX_FIELD_LEN)

    project = models.ForeignKey(Project, null=True)

    uid = models.CharField(max_length=32, unique=True)

    def save(self, *args, **kwargs):
        self.uid = self.uid or util.get_uuid(8)

        super(DataType, self).save(*args, **kwargs)

    def __str__(self):
        return self.name


class Access(models.Model):
    """
    Allows access of users to Projects.
    """

    # The numerical values for permissions matter!
    # A higher number implies all lesser permissions.
    # READ_ACCESS < EXECUTE_ACCESS < ADMIN_ACCESS
    NO_ACCESS, READ_ACCESS, RECIPE_ACCESS, EXECUTE_ACCESS, EDIT_ACCESS, UPLOAD_ACCESS, ADMIN_ACCESS = range(1, 8)
    ACCESS_CHOICES = [
        (NO_ACCESS, "No Access"),
        (READ_ACCESS, "Read Access"),
        (RECIPE_ACCESS, "Recipe Access"),
        (EXECUTE_ACCESS, "Execute Access"),
        (EDIT_ACCESS, "Edit Access"),
        (UPLOAD_ACCESS, "Upload Access"),
        (ADMIN_ACCESS, "Admin Access"),
    ]
    ACCESS_MAP = dict(ACCESS_CHOICES)

    user = models.ForeignKey(User)
    project = models.ForeignKey(Project)
    access = models.IntegerField(default=READ_ACCESS, choices=ACCESS_CHOICES)
    date = models.DateTimeField(auto_now_add=True)


@receiver(post_save, sender=Project)
def create_access(sender, instance, created, **kwargs):

    if created:
        # Creates an admin access for the user.
        access = Access.objects.create(user=instance.owner, project=instance, access=Access.ADMIN_ACCESS)
        access.save()


@receiver(post_save, sender=Project)
def add_datatypes(sender, instance, created, **kwargs):

    # Add all constant datatypes to a project on creation

    for numeric, symbol, name in const.DATA_TUPLES:

        if not DataType.objects.filter(project=instance, numeric=numeric).first():

            datatype = DataType.objects.create(project=instance, numeric=numeric,
                                    symbol=symbol, name=name)
            datatype.save()


class Data(models.Model):
    PENDING, READY, ERROR, DELETED = 1, 2, 3, 4
    STATE_CHOICES = [(PENDING, "Pending"), (READY, "Ready"), (ERROR, "Error"), (DELETED, "Deleted")]
    state = models.IntegerField(default=PENDING, choices=STATE_CHOICES)

    name = models.CharField(max_length=MAX_NAME_LEN, default="name")
    summary = models.TextField(default='summary', blank=True, max_length=MAX_TEXT_LEN)
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN)
    sticky = models.BooleanField(default=False)

    owner = models.ForeignKey(User, null=True)
    text = models.TextField(default='description', max_length=MAX_TEXT_LEN, blank=True)
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    data_type = models.IntegerField(default=GENERIC_TYPE)
    project = models.ForeignKey(Project)
    size = models.IntegerField(default=0)

    # FilePathField points to an existing file
    file = models.FilePathField(max_length=MAX_FIELD_LEN)

    uid = models.CharField(max_length=32)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.name = self.name[-MAX_NAME_LEN:]
        self.uid = self.uid or util.get_uuid(8)
        self.date = self.date or now
        self.html = make_html(self.text)
        self.owner = self.owner or self.project.owner
        self.data_type = self.data_type or GENERIC_TYPE

        # Build the data directory.
        data_dir = self.get_data_dir()
        if not os.path.isdir(data_dir):
            os.makedirs(data_dir)

        # Set the table of contents for the file.
        self.file = self.get_path()

        # Make this file if it does not exist
        if not os.path.isfile(self.file):
            with open(self.file, 'wt') as fp:
                pass

        super(Data, self).save(*args, **kwargs)

    def peek(self):
        """
        Returns a preview of the data
        """
        target = self.get_path()
        lines = open(target, 'rt').readlines()
        if len(lines) == 1:
            target = lines[0]
            return util.smart_preview(target)
        else:
            data_dir = self.get_data_dir()
            rels = [os.path.relpath(path, data_dir) for path in lines]
            return "".join(rels)

    def __str__(self):
        return self.name

    def get_data_dir(self):
        "The data directory"
        return join(self.get_project_dir(), f"store-{self.uid}")

    def get_project_dir(self):
        return self.project.get_project_dir()

    def get_path(self):
        return join(self.get_data_dir(), f"toc-{self.uid}.txt")

    def can_unpack(self):
        cond = str(self.get_path()).endswith("tar.gz")
        return cond

    def get_files(self):
        fnames = [line.strip() for line in open(self.get_path(), 'rt')]

        return fnames if len(fnames) else [""]

    def get_url(self):
        return (reverse('data_view', kwargs=dict(id=self.id)))

    def fill_dict(self, obj):
        """
        Mutates a dictionary to add more information.
        """
        fnames = self.get_files()
        if fnames:
            obj['value'] = fnames[0]
        else:
            obj['value'] = 'MISSING'

        obj['files'] = fnames
        obj['toc'] = self.get_path()
        obj['id'] = self.id
        obj['name'] = self.name
        obj['uid'] = self.uid
        obj['data_dir'] = self.get_data_dir()
        obj['project_dir'] = self.get_project_dir()
        obj['data_url'] = self.get_url()



class Analysis(models.Model):
    AUTHORIZED, UNDER_REVIEW = 1, 2

    AUTH_CHOICES = [(AUTHORIZED, "Authorized"), (UNDER_REVIEW, "Authorization Required")]

    ACTIVE, DELETED = 1, 2
    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted")]

    uid = models.CharField(max_length=32, unique=True)
    sticky = models.BooleanField(default=False)
    name = models.CharField(max_length=MAX_NAME_LEN, default="No name")
    summary = models.TextField(default='No summary.')
    text = models.TextField(default='No description.', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    owner = models.ForeignKey(User)

    security = models.IntegerField(default=UNDER_REVIEW, choices=AUTH_CHOICES)
    state = models.IntegerField(default=ACTIVE, choices=STATE_CHOICES)

    project = models.ForeignKey(Project)

    json_text = models.TextField(default="{}", max_length=MAX_TEXT_LEN)
    template = models.TextField(default="")

    date = models.DateTimeField(auto_now_add=True, blank=True)
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN, help_text="Optional image")

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
        self.name = self.name[:MAX_NAME_LEN]
        self.html = make_html(self.text)

        super(Analysis, self).save(*args, **kwargs)

    def get_project_dir(self):
        return self.project.get_project_dir()

    def url(self):
        return reverse("recipe_view", kwargs=dict(id=self.id))

    def authorized(self):
        return self.security == self.AUTHORIZED

class Job(models.Model):
    AUTHORIZED, UNDER_REVIEW = 1, 2
    AUTH_CHOICES = [(AUTHORIZED, "Authorized"), (UNDER_REVIEW, "Authorization Required")]

    QUEUED, RUNNING, COMPLETED, ERROR, DELETED, SPOOLED, PAUSED, ZOMBIE = range(1, 9)

    STATE_CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"), (PAUSED, "Paused"),
                     (SPOOLED, "Spooled"), (COMPLETED, "Completed"),
                     (ERROR, "Error"), (DELETED, "Deleted"), (ZOMBIE, "Zombie")]

    name = models.CharField(max_length=MAX_NAME_LEN, default="name")
    summary = models.TextField(default='summary')
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN)

    owner = models.ForeignKey(User)
    text = models.TextField(default='description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')

    # Job creation date
    date = models.DateTimeField(auto_now_add=True)

    # Job runtime date.
    start_date = models.DateTimeField(null=True, blank=True)
    end_date = models.DateTimeField(null=True, blank=True)

    sticky = models.BooleanField(default=False)
    analysis = models.ForeignKey(Analysis)
    project = models.ForeignKey(Project)
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

    state = models.IntegerField(default=1, choices=STATE_CHOICES)

    path = models.FilePathField(default="")

    def is_running(self):
        return self.state == Job.RUNNING

    def __str__(self):
        return self.name

    def get_url(self, path=''):
        "Return the url to the job directory"
        return f"jobs/job-{self.uid}/" + path

    def url(self):
        return reverse("job_view", kwargs=dict(id=self.id))

    def get_project_dir(self):
        return self.project.get_project_dir()

    @property
    def json_data(self):
        "Returns the json_text as parsed json_data"
        return hjson.loads(self.json_text)

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

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.name = self.name or self.analysis.name
        self.date = self.date or now
        self.html = make_html(self.text)
        self.name = self.name[:MAX_NAME_LEN]
        self.uid = self.uid or util.get_uuid(8)
        self.template = self.analysis.template
        self.security = self.analysis.security
        self.stderr_log = self.stderr_log[:MAX_LOG_LEN]
        self.stdout_log = self.stdout_log[:MAX_LOG_LEN]
        self.name = self.name or self.analysis.name
        # write an index.html to the file
        if not os.path.isdir(self.path):
            path = join(settings.MEDIA_ROOT, "jobs", f"job-{self.uid}")
            os.makedirs(path)
            self.path = path

        super(Job, self).save(*args, **kwargs)
