import mimetypes
from itertools import islice

import hjson as json
import mistune
from django.contrib.auth.models import User
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.template.loader import get_template

from django.contrib.auth.models import Group
import mistune

from django.urls import reverse
from django.utils import timezone
from . import settings
from . import util
from .const import *

def join(*args):
    return os.path.abspath(os.path.join(*args))


def make_html(text):
    return mistune.markdown(text)


def get_datatype(file):
    return Data.FILE


def directory_path(instance, filename):
    pieces = os.path.basename(filename).split(".")
    # May have multiple extensions
    exts = ".".join(pieces[1:]) or "data"
    uid = util.get_uuid(8)
    filename = f"data-{uid}.{exts}"

    return f'{instance.project.get_path()}/{filename}'


class Project(models.Model):

    title = models.CharField(max_length=256)
    owner = models.ForeignKey(User)
    text = models.TextField(default='text')
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    # Project restircted to one group
    group = models.ForeignKey(Group)

    uid = models.CharField(max_length=32)

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

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
        return self.title

    def url(self):
        return reverse("project_view", kwargs=dict(id=self.id))

    def get_path(self):
        return join(settings.MEDIA_ROOT, "projects", f"proj-{self.uid}")

    def create_analysis(self, json_text, template, owner=None, summary='', title='', text=''):
        """
        Creates analysis from a spec and template
        """
        owner = owner or self.owner
        title = title or 'Analysis Title'
        text = text or 'Analysis Text'
        analysis = Analysis.objects.create(project=self, summary=summary, json_text=json_text,
                                           owner=owner, title=title, text=text,
                                           template=template, )
        return analysis


class Data(models.Model):
    FILE, COLLECTION = 1, 2
    TYPE_CHOICES = [(FILE, "File"), (COLLECTION, "Collection")]

    title = models.CharField(max_length=256)
    owner = models.ForeignKey(User)
    text = models.TextField(default='text')
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)
    type = models.IntegerField(default=FILE, choices=TYPE_CHOICES)
    data_type = models.IntegerField(default=GENERIC_TYPE)
    project = models.ForeignKey(Project)
    size = models.CharField(null=True, max_length=256)

    file = models.FileField(null=True, upload_to=directory_path)
    path = models.FilePathField(null=True)

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.date = self.date or now
        self.html = make_html(self.text)
        super(Data, self).save(*args, **kwargs)

    def peek(self):
        """Peeks at the data if it is text"""
        mimetype, mimecode = mimetypes.guess_type(self.get_path())
        if mimetype == 'text/plain':
            stream = open(self.file.path)
            lines = [line for line in islice(stream, 10)]
            content = "\n".join(lines)
            return content

        return "*** Binary file ***"

    def __str__(self):
        return self.title

    def get_path(self):
        return self.file.path


class Analysis(models.Model):

    project = models.ForeignKey(Project)
    title = models.CharField(max_length=256)
    owner = models.ForeignKey(User)
    text = models.TextField(default='text')
    summary = models.TextField(default='summary')
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)
    json_text = models.TextField(default="{}")
    template = models.TextField(default="makefile")

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

    def __str__(self):
        return self.title

    @property
    def json_data(self):
        "Returns the json_text as parsed json_data"
        return json.loads(self.json_text)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.date = self.date or now
        self.html = make_html(self.text)
        super(Analysis, self).save(*args, **kwargs)

    def create_job(self, json_text='', json_data={}, owner=None, title=None, state=None):
        """
        Creates a job from an analysis.
        """
        title = title or self.title
        state = state or Job.QUEUED
        owner = owner or self.project.owner

        if json_data:
            json_text = json.dumps(json_data)
        else:
            json_text = json_text or self.json_text

        job = Job.objects.create(title=title, summary=self.summary, state=state, json_text=json_text,
                                 project=self.project, analysis=self, owner=owner,
                                 template=self.template)

        return job


class Job(models.Model):
    # file path to media
    QUEUED, RUNNING, FINISHED, ERROR = 1, 2, 3, 4

    STATE_CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"),
                     (FINISHED, "Finished"), (ERROR, "Error")]


    title = models.CharField(max_length=256)
    owner = models.ForeignKey(User)
    text = models.TextField(default='text')
    summary = models.TextField(default='summary')
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    analysis = models.ForeignKey(Analysis)
    project = models.ForeignKey(Project)
    json_text = models.TextField(default="commands")

    uid = models.CharField(max_length=32)
    template = models.TextField(default="makefile")
    log = models.TextField(default="No data logged for current job")

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

    state = models.IntegerField(default=1, choices=STATE_CHOICES)

    path = models.FilePathField(default="")

    def is_running(self):
        return self.state == Job.RUNNING

    def __str__(self):
        return self.title

    def get_url(self, path=''):
        "Return the url to the job directory"
        return f"jobs/job-{self.uid}/" + path

    @property
    def json_data(self):
        "Returns the json_text as parsed json_data"
        return json.loads(self.json_text)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.date = self.date or now
        self.html = make_html(self.text)

        self.uid = self.uid or util.get_uuid(8)
        self.template = self.analysis.template

        self.title = self.title or self.analysis.title
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
        instance.groups.add(Group.objects.get(name='Public'))


post_save.connect(create_profile, sender=User)
