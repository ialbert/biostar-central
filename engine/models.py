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

    project = models.ForeignKey(Project)
    #pipeline = ""

    def save(self, *args, **kwargs):
        super(Analysis, self).save(*args, **kwargs)


class Result(Base):

    QUEUED, RUNNING, FINISHED, ERROR = 1,2,3,4
    CHOICES = [("Queued", QUEUED), ("Running", RUNNING),
               ("Finished", FINISHED), ("Error", ERROR)]
    analysis = models.ForeignKey(Analysis)
    state =  models.IntegerField(default=1, choices=CHOICES)
    directory = models.FilePathField(default="media")
    commands = models.TextField(default="commands")

    def save(self, *args, **kwargs):
        super(Result, self).save(*args, **kwargs)

    
class Profile(models.Model):

    user = models.ForeignKey(User)


@receiver(post_save, sender=User)
def create_profile(sender, instance, created, **kwargs):
    if created:
        Profile.objects.create(user=instance)


#post_save.connect(create_profile, sender=User)
