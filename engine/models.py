import logging, hjson

from django.contrib.auth.models import User
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver


from django.contrib.auth.models import Group
import mistune

from django.urls import reverse
from django.utils import timezone
from . import settings
from . import util
from .const import *
from biostar.tools import defaults
from django.core.files import File
logger = logging.getLogger("engine")

# The maximum length in characters for a typical name and text field.
MAX_NAME_LEN = 256
MAX_TEXT_LEN = 10000

def join(*args):
    return os.path.abspath(os.path.join(*args))


def make_html(text):
    return mistune.markdown(text)


def get_datatype(file):
    return Data.FILE


def filter_by_usage():
    return


def upload_path(instance, filename):
    # Name the data by the filename.
    pieces = os.path.basename(filename).split(".")
    # File may have multiple extensions
    exts = ".".join(pieces[1:]) or "data"
    dataname = f"data-{instance.uid}.{exts}"

    return join(instance.project.get_path(), f"data-{instance.uid}", dataname)


class Project(models.Model):

    ADMIN, USER = 1, 2
    USAGE_CHOICES = [(ADMIN, "admin"), (USER, "user")]
    usage = models.IntegerField(default=USER, choices=USAGE_CHOICES)

    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')

    owner = models.ForeignKey(User)
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)

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
        return self.name

    def url(self):
        return reverse("project_view", kwargs=dict(id=self.id))

    def get_path(self):
        return join(settings.MEDIA_ROOT, "projects", f"proj-{self.uid}")

    def get_data(self, data_type=None):
        """
        Returns a dictionary keyed by data stored in the project.
        """
        query = Data.objects.filter(project=self)
        if data_type:
            query = query.filter(data_type=data_type)
        datamap = dict( (obj.id, obj) for obj in query )
        return datamap

    def create_analysis(self, json_text, template, owner=None, summary='', name='', text='', usage=None):
        """
        Creates analysis from a spec and template
        """
        owner = owner or self.owner
        name = name or 'Analysis name'
        text = text or 'Analysis text'
        usage = usage or defaults.USAGE
        analysis = Analysis.objects.create(project=self, summary=summary, json_text=json_text,
                                           owner=owner, name=name, text=text, usage=usage,
                                           template=template )

        logger.info(f"Created analysis id={analysis.id} of usage_type={dict(Analysis.USAGE_CHOICES)[analysis.usage]}")
        return analysis

    def create_data(self,  stream=None, fname=None, name="data.bin", owner=None, text='', data_type=None, usage=None):
        """
        Creates a data for the project from filename or a stream.
        """

        if fname:
            stream = File(open(fname, 'rb'))
            name = os.path.basename(fname)
        owner = owner or self.owner
        text = text or "No description"
        data_type = data_type or GENERIC_TYPE
        usage = usage or defaults.USAGE
        data = Data(name=name, owner=owner, usage=usage,
                    text=text, project=self, data_type=data_type)

        # Need to save before uid gets triggered.
        data.save()
        # This saves the into the
        data.file.save(name, stream, save=True)

        # Set the pending to ready after the file saves.
        Data.objects.filter(id=data.id).update(state=Data.READY)

        # Updates its own size.
        data.set_size()

        logger.info(f"Added data id={data.id} of type={data.data_type}")
        return data


class Data(models.Model):

    ADMIN, USER = 1, 2
    FILE, COLLECTION = 1, 2
    PENDING, READY = 1,2

    TYPE_CHOICES = [(FILE, "File"), (COLLECTION, "Collection")]
    USAGE_CHOICES = [(ADMIN, "Admin"), (USER, "User")]
    STATE_CHOICES = [(PENDING, "Pending"), (READY, "Ready")]

    usage = models.IntegerField(default=USER, choices=USAGE_CHOICES)
    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')

    owner = models.ForeignKey(User)
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)
    type = models.IntegerField(default=FILE, choices=TYPE_CHOICES)

    data_type = models.IntegerField(default=GENERIC_TYPE)
    project = models.ForeignKey(Project)
    size = models.CharField(null=True, max_length=256)

    state = models.IntegerField(default=PENDING, choices=STATE_CHOICES)
    file = models.FileField(null=True, upload_to=upload_path, max_length=500)
    uid = models.CharField(max_length=32)

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.uid = self.uid or util.get_uuid(8)
        self.date = self.date or now
        self.html = make_html(self.text)

        if not os.path.isdir(join(self.project.get_path(), f"data-{self.uid}")):
            os.makedirs(join(self.project.get_path(), f"data-{self.uid}"))

        super(Data, self).save(*args, **kwargs)

    def peek(self):
        """
        Returns a preview of the data
        """
        return util.smart_preview(self.get_path())

    def set_size(self):
        """
        Sets the size of the data.
        """
        try:
            size = os.path.getsize(self.get_path())
        except:
            size = 0
        Data.objects.filter(id=self.id).update(size=size)
    def __str__(self):
        return self.name

    def get_path(self):
        return self.file.path

    def fill_dict(self, obj):
        """
        Mutates a dictionary to add more information.
        """
        obj['path'] = self.get_path()
        obj['value'] = self.id
        obj['name'] = self.name


class Analysis(models.Model):

    ADMIN, USER = 1, 2
    USAGE_CHOICES = [(ADMIN, "admin"), (USER, "user")]
    usage = models.IntegerField(default=USER, choices=USAGE_CHOICES)

    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    owner = models.ForeignKey(User)

    project = models.ForeignKey(Project)

    json_text = models.TextField(default="{}")
    template = models.TextField(default="makefile")
    uid = models.CharField(max_length=32)

    date = models.DateTimeField(auto_now_add=True, blank=True)

    # Job start and end.
    start_date = models.DateTimeField(null=True, blank=True)
    end_date = models.DateTimeField(null=True, blank=True)

    # Will be false if the object is deleted.
    valid = models.BooleanField(default=True)

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

        self.html = make_html(self.text)
        super(Analysis, self).save(*args, **kwargs)

    def create_job(self, json_text='', json_data={}, owner=None, name=None, state=None, usage=None):
        """
        Creates a job from an analysis.
        """
        name = name or self.name
        state = state or Job.QUEUED
        owner = owner or self.project.owner
        usage = usage or defaults.USAGE

        if json_data:
            json_text = hjson.dumps(json_data)
        else:
            json_text = json_text or self.json_text

        job = Job.objects.create(name=name, summary=self.summary, state=state, json_text=json_text,
                                 project=self.project, analysis=self, owner=owner, usage=usage,
                                 template=self.template)

        logger.info(f"Queued job: '{job.name}'")
        return job


class Job(models.Model):

    ADMIN, USER = 1, 2
    QUEUED, RUNNING, FINISHED, ERROR = 1, 2, 3, 4
    USAGE_CHOICES = [(ADMIN, "admin"), (USER, "user")]
    STATE_CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"),
                     (FINISHED, "Finished"), (ERROR, "Error")]

    usage = models.IntegerField(default=USER, choices=USAGE_CHOICES)
    name = models.CharField(max_length=256, default="no name")
    summary = models.TextField(default='no summary')

    owner = models.ForeignKey(User)
    text = models.TextField(default='no description', max_length=MAX_TEXT_LEN)
    html = models.TextField(default='html')
    date = models.DateTimeField(auto_now_add=True)

    analysis = models.ForeignKey(Analysis)
    project = models.ForeignKey(Project)
    json_text = models.TextField(default="commands")

    uid = models.CharField(max_length=32)
    template = models.TextField(default="makefile")

    # This will be set when the job attempts to run.
    script = models.TextField(default="")

    # Keeps track of errors.
    log = models.TextField(default="No data logged for current job", max_length=10 * MAX_TEXT_LEN)

    # Will be false if the objects is to be deleted.
    valid = models.BooleanField(default=True)

    state = models.IntegerField(default=1, choices=STATE_CHOICES)

    path = models.FilePathField(default="")
    local = models.FilePathField(default="")

    def is_running(self):
        return self.state == Job.RUNNING

    def __str__(self):
        return self.name

    def get_url(self, path=''):
        "Return the url to the job directory"
        return f"jobs/job-{self.uid}/" + path

    @property
    def json_data(self):
        "Returns the json_text as parsed json_data"
        return hjson.loads(self.json_text)

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.date = self.date or now
        self.html = make_html(self.text)

        self.uid = self.uid or util.get_uuid(8)
        self.template = self.analysis.template

        self.name = self.name or self.analysis.name
        # write an index.html to the file
        if not os.path.isdir(self.path):

            path = join(settings.MEDIA_ROOT, "jobs", f"job-{self.uid}")
            local = join(settings.LOCAL_ROOT, "jobs", f"job-{self.uid}")

            os.makedirs(path)
            os.makedirs(local)

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
        # instance.groups.add(Group.objects.get(name='Public'))


post_save.connect(create_profile, sender=User)
