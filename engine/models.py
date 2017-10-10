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


class Data(Base):
    FILE, COLLECTION = 1, 2
    TYPE_CHOICES =[(FILE, "File"),(COLLECTION ,"Collection")]
    type = models.IntegerField(default=FILE, choices=TYPE_CHOICES)

    project = models.ForeignKey(Project)
    file = models.FileField(null=True)
    path = models.FilePathField(null=True)


    def save(self, *args, **kwargs):
        # chack if file or path is set or raise error
        super(Data, self).save(*args, **kwargs)


class Analysis(Base):

    spec_origin = models.TextField(default="media")
    spec_source = models.TextField(default=r"""{ 
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

    def save(self, *args, **kwargs):

        # Save source from a string if one is given.
        if kwargs.get("spec_source"):
            source = kwargs.pop("spec_source")
            self.spec_source = source
        else:
            self.spec_source = json.dumps(util.safe_load(self.spec_origin))

        super(Analysis, self).save(*args, **kwargs)


class Job(Base):

    # file path to media
    QUEUED, RUNNING, FINISHED, ERROR = 1, 2, 3, 4
    STATE_CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"),
               (FINISHED, "Finished"), (ERROR, "Stopped")]
    STATE_MAP = dict(STATE_CHOICES)

    analysis = models.ForeignKey(Analysis)
    project = models.ForeignKey(Project)
    json_data = models.TextField(default="commands")

    # uniqe directory creation
    uid = models.CharField(max_length=32)
    makefile_template = models.TextField(default="makefile")
    log = models.TextField(default="log")

    state = models.IntegerField(default=1, choices=STATE_CHOICES)
    # rename to path
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
