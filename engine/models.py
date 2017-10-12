import hjson as json
import os
from django.db import models
from django.contrib.auth.models import User
from django.utils import timezone
from django.db.models.signals import post_save
from django.dispatch import receiver
import mistune
from . import settings
from . import util
from django.urls import reverse


def join(*args):
    return os.path.abspath(os.path.join(*args))


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
    project = models.ForeignKey(Project)
    size = ''

    file = models.FileField(null=True, upload_to=directory_path)
    path = models.FilePathField(null=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def save(self, *args, **kwargs):
        super(Data, self).save(*args, **kwargs)

    def get_path(self):
        return self.file if self.type == Data.FILE else self.path

class Analysis(Base):

    json_file = models.TextField(default="media")
    json_data = models.TextField(default=r"""{ 
                                            field1_name: 
                                            {
                                            value:field_value, 
                                            type:FIELD_TYPE
                                            }
                                            field2_name:
                                            {
                                            value:field_value, 
                                            type:FIELD_TYPE
                                            }
                                            }""")

    makefile_template = models.TextField(default="makefile")

    def save(self, *args, **kwargs):
        super(Analysis, self).save(*args, **kwargs)


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
    # jobstart and job_end fields.
    state = models.IntegerField(default=1, choices=STATE_CHOICES)
    path = models.FilePathField(default="")

    def save(self, *args, **kwargs):

        self.uid = self.uid or util.get_uuid(8)
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
