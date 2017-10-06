
import os
from django.db import models
from django.contrib.auth.models import User
from django.utils import timezone
from django.db.models.signals import post_save
from django.dispatch import receiver
import mistune
from .util import get_uuid
from django.urls import reverse

def make_html(text):

    return mistune.markdown(text)


class Base(models.Model):

    title = models.CharField(max_length=256)
    owner = models.ForeignKey(User)
    text = models.TextField(default='text')
    html = models.TextField(default='html')
    date = models.DateTimeField()

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

    def save(self, *args, **kwargs):
        super(Project, self).save(*args, **kwargs)

    def url(self):
        return reverse("project_view", kwargs=dict(id=self.id))


class Data(Base):

    project = models.ForeignKey(Project)
    file = models.FileField(null=True)

    def save(self, *args, **kwargs):
        super(Data, self).save(*args, **kwargs)


class Analysis(Base):

    #project = models.ForeignKey(Project)
    json_spec = models.TextField(default="# Enter json commands")
    # = models.FileField(default="media")
    #makefile_template = models.TextField(default="media")

    def save(self, *args, **kwargs):
        super(Analysis, self).save(*args, **kwargs)

    # set json_spec from json_file
    # assumes json_file is already a json file

class Job(Base):

    # file path to media
    QUEUED, RUNNING, FINISHED, ERROR = 1, 2, 3, 4
    CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"),
               (FINISHED, "Finished"), (ERROR, "Error")]

    analysis = models.ForeignKey(Analysis)
    project = models.ForeignKey(Project)
    json_data = models.TextField(default="commands")
    uid = models.CharField(default="", max_length=32)
    #makefile_template = models.TextField(default="media")
    makefile_template = models.TextField(default="makefile")
    log = models.TextField(default="log")

    state = models.IntegerField(default=1, choices=CHOICES)
    directory = models.FilePathField(default="media")


    def save(self, *args, **kwargs):
        # create
        self.uid = self.uid or get_uuid()
        if not os.path.isdir(self.directory):
            path = os.path.abspath(os.path.join(settings.MEDIA_DIR, uuid))
            os.mkdir(path)
            self.directory = path
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
