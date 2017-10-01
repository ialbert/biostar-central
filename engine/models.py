import json

from django.db import models
from django.contrib.auth.models import User
from django.utils import timezone
from django.db.models.signals import post_save
from django.dispatch import receiver
import mistune


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


class Data(Base):

    project = models.ForeignKey(Project)
    file = models.FileField(null=True)

    def save(self, *args, **kwargs):
        super(Data, self).save(*args, **kwargs)


class Analysis(Base):

    #project = models.ForeignKey(Project)
    json_spec = models.TextField(default="commands")
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

    #makefile_template = models.TextField(default="media")
    makefile = models.TextField(default="media")

    state = models.IntegerField(default=1, choices=CHOICES)

    def save(self, *args, **kwargs):
        super(Job, self).save(*args, **kwargs)

    def change_state(self):
        return NotImplemented



class Result(Base):

    # file path to media

    job = models.ForeignKey(Job)
    directory = models.FilePathField(default="media")

    def save(self, *args, **kwargs):
        super(Result, self).save(*args, **kwargs)


class Profile(models.Model):

    user = models.ForeignKey(User)


@receiver(post_save, sender=User)
def create_profile(sender, instance, created, **kwargs):
    if created:
        Profile.objects.create(user=instance)


#post_save.connect(create_profile, sender=User)
