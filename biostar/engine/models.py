import logging
import hjson
import mistune
from django.db import models
from django.db.models.signals import post_save, pre_save
from django.dispatch import receiver
from django.urls import reverse
from django.utils import timezone

from biostar import settings
from biostar.accounts.models import User
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

    # Affects the sort order.
    sticky = models.BooleanField(default=False)

    # Limits who can access the project.
    privacy = models.IntegerField(default=PRIVATE, choices=PRIVACY_CHOICES)

    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN)
    name = models.CharField(default="name", max_length=MAX_NAME_LEN)
    summary = models.TextField(default='summary', max_length=MAX_TEXT_LEN)
    deleted = models.BooleanField(default=False)

    # We need to keep the owner.
    owner = models.ForeignKey(User, null=False, on_delete=models.CASCADE)
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

    def uid_is_set(self):
        assert bool(self.uid.strip()), "Sanity check. UID should always be set."

    def url(self):
        self.uid_is_set()
        return reverse("project_view", kwargs=dict(uid=self.uid))

    def get_project_dir(self):
        self.uid_is_set()
        return join(settings.MEDIA_ROOT, "projects", f"proj-{self.uid}")


    def get_data_dir(self):
        "Match consistency of data dir calls"
        return self.get_project_dir()

    @property
    def project(self):
        return self

class Access(models.Model):
    """
    Allows access of users to Projects.
    """

    # The numerical values for permissions matter!
    # A higher number implies all lesser permissions.
    # READ_ACCESS < WRITE_ACCESS < OWNER_ACCESS
    NO_ACCESS, READ_ACCESS, WRITE_ACCESS, OWNER_ACCESS = range(1, 5)
    ACCESS_CHOICES = [
        (NO_ACCESS, "No Access"),
        (READ_ACCESS, "Read Access"),
        (WRITE_ACCESS, "Write Access"),
        (OWNER_ACCESS, "Owner Access")
    ]

    ACCESS_MAP = dict(ACCESS_CHOICES)

    user = models.ForeignKey(User, on_delete=models.CASCADE)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    access = models.IntegerField(default=READ_ACCESS, choices=ACCESS_CHOICES)
    date = models.DateTimeField(auto_now_add=True)

    def save(self, *args, **kwargs):
        super(Access, self).save(*args, **kwargs)



@receiver(post_save, sender=Project)
def update_access(sender, instance, created, raw, update_fields, **kwargs):

    # Drop previous OWNER_ACCES permissions if these exists.
    # This is needed when projects change owners.
    Access.objects.filter(project=instance, access=Access.OWNER_ACCESS).delete()

    # Create the admin access for the current owner.
    access = Access.objects.create(user=instance.owner, project=instance, access=Access.OWNER_ACCESS)
    access.save()

    # Checking the name is implemented as a post-save signal to
    i = 0
    name = instance.name
    # Count works better since .exist() will include the first upload
    while Project.objects.filter(name=name).count() > 1 and i < 100:
        i += 1
        name = f"{instance.name} ({i})"
        Project.objects.filter(uid=instance.uid, pk=instance.pk).update(name=name)



class Data(models.Model):
    PENDING, READY, ERROR, = 1, 2, 3
    STATE_CHOICES = [(PENDING, "Pending"), (READY, "Ready"), (ERROR, "Error")]
    state = models.IntegerField(default=PENDING, choices=STATE_CHOICES)

    LINK, UPLOAD = 1, 2
    METHOD_CHOICE = [(LINK, "Linked Data"), (UPLOAD, "Uploaded Data")]
    method = models.IntegerField(default=LINK, choices=METHOD_CHOICE)

    name = models.CharField(max_length=MAX_NAME_LEN, default="name")
    summary = models.TextField(default='summary', blank=True, max_length=MAX_TEXT_LEN)
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN)
    sticky = models.BooleanField(default=False)
    deleted = models.BooleanField(default=False)

    owner = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    text = models.TextField(default='description', max_length=MAX_TEXT_LEN, blank=True)
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    type = models.CharField(max_length=MAX_NAME_LEN, default="DATA")
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    size = models.IntegerField(default=0)

    # FilePathField points to an existing file
    file = models.FilePathField(max_length=MAX_FIELD_LEN)

    uid = models.CharField(max_length=32, unique=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.name = self.name[-MAX_NAME_LEN:]
        self.uid = self.uid or util.get_uuid(8)
        self.date = self.date or now
        self.html = make_html(self.text)
        self.owner = self.owner or self.project.owner


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
        assert self.uid, "Sanity check. UID should always be set."
        return join(self.get_project_dir(), f"store-{self.uid}")

    def get_project_dir(self):
        return self.project.get_project_dir()

    def get_path(self):
        path = join(settings.TOC_ROOT, f"toc-{self.uid}.txt")
        return path

    def can_unpack(self):
        cond = str(self.get_path()).endswith("tar.gz")
        return cond

    def get_files(self):
        fnames = [line.strip() for line in open(self.get_path(), 'rt')]
        return fnames if len(fnames) else [""]

    def get_url(self, path=""):
        "Returns url to the data directory"
        return f"projects/proj-{self.project.uid}/store-{self.uid}/" + path

    def url(self):
        return (reverse('data_view', kwargs=dict(uid=self.uid)))

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
        obj['data_url'] = self.url()


@receiver(post_save, sender=Data)
def check_data_name(sender, instance, created, **kwargs):
    # Add primary key to the name if it already exists.

    i = 0
    name = instance.name
    # Count works better since .exists() includes first uploads
    # causing everything to start with a (1). 
    while Data.objects.filter(name=name, project=instance.project).count() > 1 and i < 100:
        i += 1
        name = f"{instance.name} ({i})"
        Data.objects.filter(project=instance.project, uid=instance.uid, pk=instance.pk).update(name=name)


class Analysis(models.Model):
    AUTHORIZED, UNDER_REVIEW = 1, 2
    AUTH_CHOICES = [(AUTHORIZED, "Authorized"), (UNDER_REVIEW, "Authorization Required")]
    security = models.IntegerField(default=UNDER_REVIEW, choices=AUTH_CHOICES)

    deleted = models.BooleanField(default=False)
    uid = models.CharField(max_length=32, unique=True)
    sticky = models.BooleanField(default=False)
    name = models.CharField(max_length=MAX_NAME_LEN, default="No name")
    summary = models.TextField(default='No summary.')
    text = models.TextField(default='No description.', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    owner = models.ForeignKey(User, on_delete=models.CASCADE)

    diff_author = models.ForeignKey(User, on_delete=models.CASCADE, related_name="diff_author", null=True)
    diff_date = models.DateField(blank=True, auto_now_add=True)

    project = models.ForeignKey(Project, on_delete=models.CASCADE)

    json_text = models.TextField(default="{}", max_length=MAX_TEXT_LEN)
    template = models.TextField(default="")
    last_valid = models.TextField(default='')

    date = models.DateTimeField(auto_now_add=True, blank=True)
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN,
                              help_text="Optional image")

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
        self.diff_date = self.diff_date or now
        self.name = self.name[:MAX_NAME_LEN]
        self.html = make_html(self.text)
        self.diff_author = self.diff_author or self.owner

        # Ensure Unix line endings.
        self.template = self.template.replace('\r\n', '\n') if self.template else ""

        if self.security == self.AUTHORIZED:
            self.last_valid = self.template

        super(Analysis, self).save(*args, **kwargs)

    def get_project_dir(self):
        return self.project.get_project_dir()

    def url(self):
        assert self.uid, "Sanity check. UID should always be set."
        return reverse("recipe_view", kwargs=dict(uid=self.uid))

    def authorized(self):
        return self.security == self.AUTHORIZED


class Job(models.Model):
    AUTHORIZED, UNDER_REVIEW = 1, 2
    AUTH_CHOICES = [(AUTHORIZED, "Authorized"), (UNDER_REVIEW, "Authorization Required")]

    QUEUED, RUNNING, COMPLETED, ERROR, SPOOLED, PAUSED = range(1, 7)

    STATE_CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"), (PAUSED, "Paused"),
                     (SPOOLED, "Spooled"), (COMPLETED, "Completed"), (ERROR, "Error")]

    state = models.IntegerField(default=QUEUED, choices=STATE_CHOICES)

    deleted = models.BooleanField(default=False)
    name = models.CharField(max_length=MAX_NAME_LEN, default="name")
    summary = models.TextField(default='summary')
    image = models.ImageField(default=None, blank=True, upload_to=image_path, max_length=MAX_FIELD_LEN)

    owner = models.ForeignKey(User, on_delete=models.CASCADE)
    text = models.TextField(default='description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')

    # Job creation date
    date = models.DateTimeField(auto_now_add=True)

    # Job runtime date.
    start_date = models.DateTimeField(null=True, blank=True)
    end_date = models.DateTimeField(null=True, blank=True)

    sticky = models.BooleanField(default=False)
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

    def is_running(self):
        return self.state == Job.RUNNING

    def __str__(self):
        return self.name

    def get_url(self, path=''):
        "Return the url to the job directory"

        return f"jobs/job-{self.uid}/" + path

    def url(self):
        return reverse("job_view", kwargs=dict(uid=self.uid))

    def get_project_dir(self):
        return self.project.get_project_dir()

    def get_data_dir(self):
        return self.path

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
        self.name = self.name or f"Results for: {self.analysis.name}"
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
