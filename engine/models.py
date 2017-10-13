import hjson as json
import os
from django.db import models
from django.contrib.auth.models import User
from django.utils import timezone
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.template.loader import get_template
import mistune
from . import settings
from . import util
from django.urls import reverse
from .const import *

def join(*args):
    return os.path.abspath(os.path.join(*args))


def make_html(text):
    return mistune.markdown(text)


class Base(models.Model):

    title = models.CharField(max_length=256)
    owner = models.ForeignKey(User)
    text = models.TextField(default='text')
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.title

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.date = self.date or now
        self.html = make_html(self.text)
        super(Base, self).save(*args, **kwargs)

    class Meta:
        abstract = True


class Project(Base):

    uid = models.CharField(max_length=32)

    def save(self, *args, **kwargs):
        self.uid = self.uid or util.get_uuid(8)
        if not os.path.isdir(self.get_path()):
            os.makedirs(self.get_path(), exist_ok=True)

        super(Project, self).save(*args, **kwargs)

    def url(self):
        return reverse("project_view", kwargs=dict(id=self.id))

    def get_path(self):
        return join(settings.DATA_ROOT, self.uid)

def directory_path(instance, filename):

    return f'{instance.project.get_path()}/{filename}'


class Data(Base):

    FILE, COLLECTION = 1, 2
    TYPE_CHOICES =[(FILE, "File"),(COLLECTION ,"Collection")]
    type = models.IntegerField(default=FILE, choices=TYPE_CHOICES)
    data_type = models.IntegerField(default=GENERIC_TYPE)
    project = models.ForeignKey(Project)

    size = models.CharField(null=True, max_length=256)

    file = models.FileField(null=True, upload_to=directory_path)
    path = models.FilePathField(null=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def save(self, *args, **kwargs):
        #self.size = f"{models.FileField.__sizeof__(self.file)}"
        super(Data, self).save(*args, **kwargs)

    def get_path(self):
        return self.file if self.type == Data.FILE else self.path


def make_analysis_from_spec(path, user, project):

    json_obj = util.safe_load(path)
    title = json_obj["analysis_spec"]["title"]
    text = json_obj["analysis_spec"]["text"]
    template_path = json_obj["template"]["path"]

    makefile_template = get_template(template_path).template.source

    analysis = Analysis(json_data=json.dumps(json_obj), owner=user, title=title, text=text,
                        makefile_template=makefile_template, project=project)
    analysis.save()

    return analysis



class Analysis(Base):

    json_data = models.TextField(default="{}")
    makefile_template = models.TextField(default="makefile")
    project = models.ForeignKey(Project)

    def save(self, *args, **kwargs):
        super(Analysis, self).save(*args, **kwargs)


def make_job(owner, analysis, project, json_data=None, title=None, state=None):

    title = title or analysis.title
    state = state or Job.QUEUED
    filled_json = json_data or analysis.json_data

    job = Job(title=title, state=state, json_data=filled_json,
              project=project, analysis=analysis, owner=owner,
              makefile_template=analysis.makefile_template)
    job.save()

    return job


class Job(Base):

    # file path to media
    QUEUED, RUNNING, FINISHED, ERROR = 1, 2, 3, 4
    STATE_CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"),
               (FINISHED, "Finished"), (ERROR, "Stopped")]

    analysis = models.ForeignKey(Analysis)
    project = models.ForeignKey(Project)
    json_data = models.TextField(default="commands")

    uid = models.CharField(max_length=32)
    makefile_template = models.TextField(default="makefile")

    log = models.TextField(default="log")

    #job_start = models.DateField(auto_now_add=True)
    #job_killed = models.DateField()
    state = models.IntegerField(default=1, choices=STATE_CHOICES)
    path = models.FilePathField(default="")

    def save(self, *args, **kwargs):

        self.uid = self.uid or util.get_uuid(8)
        self.makefile_template = self.analysis.makefile_template

        self.title = self.title or self.analysis.title
        # write an index.html to the file
        if not os.path.isdir(self.path):
            path = os.path.abspath(os.path.join(settings.MEDIA_ROOT, self.uid))
            os.mkdir(path)
            self.path = path

        super(Job, self).save(*args, **kwargs)

    def url(self):
        return reverse("job_view", kwargs=dict(id=self.id))


class Profile(models.Model):

    user = models.ForeignKey(User)


@receiver(post_save, sender=User)
def create_profile(sender, instance, created, **kwargs):
    if created:
        Profile.objects.create(user=instance)


#post_save.connect(create_profile, sender=User)
