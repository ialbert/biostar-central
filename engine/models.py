from django.db import models
import django_tables2 as tables
from django.contrib.auth.models import User
from django.utils import timezone
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
        #self.html = make_html(self.text)
        super(Base, self).save(*args, **kwargs)

    class Meta:
        abstract = True


class Project(Base):

    def save(self, *args, **kwargs):
        super(Project, self).save(*args, **kwargs)


class Data(Base):

    project = models.ForeignKey(Project)
    #upload = ""
    # what type of file
    #filetypes = ""

    # number of files
    #nfiles = 2
    #
    #for x in range(nfiles):
    #   make a file feild

    def save(self, *args, **kwargs):
        super(Data, self).save(*args, **kwargs)


class Analysis(Base):

    project = models.ForeignKey(Project)

    def save(self, *args, **kwargs):
        super(Analysis, self).save(*args, **kwargs)


class Pipeline(Base):

    analysis = models.ForeignKey(Analysis)

    #status = ""
    def save(self, *args, **kwargs):
        super(Pipeline, self).save(*args, **kwargs)


class Result(Base):

    pipeline = models.ForeignKey(Pipeline)


    def save(self, *args, **kwargs):
        super(Result, self).save(*args, **kwargs)

    
class PipelineTable(tables.Table):

    class Meta:
        model = Pipeline
        #attrs = {"class": "table-striped table-bordered"}
        #project = Pipeline.objects.all().analysis.project.title
        fields = ("title", "date", "text",
                   "analysis", "project")

