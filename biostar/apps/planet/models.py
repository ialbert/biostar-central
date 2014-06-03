from django.db import models
from django.conf import settings
import os, urllib, logging, feedparser, datetime
from django.core.urlresolvers import reverse
from django.utils.timezone import utc
from django.contrib import admin

logger = logging.getLogger(__name__)

def now():
    return datetime.datetime.utcnow().replace(tzinfo=utc)

def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))

# Create your models here.

class Blog(models.Model):
    "Represents a blog"
    title = models.CharField(verbose_name='Blog Name', max_length=255, default="", blank=False)
    desc = models.TextField(default='', blank=True)
    feed = models.URLField()
    link = models.URLField()
    active = models.BooleanField(default=True)
    list_order = models.IntegerField(default=0)

    @property
    def fname(self):
        fname = abspath(settings.PLANET_DIR, '%s.xml' % self.id)
        return fname

    def parse(self):
        try:
            doc = feedparser.parse(self.fname)
        except Exception, exc:
            logger.error("error %s parsing blog %s", (exc, self.id))
            doc = None
        return doc

    def download(self):
        try:
            text = urllib.urlopen(self.feed).read()
            stream = file(self.fname, 'wt')
            stream.write(text)
            stream.close()
        except Exception, exc:
            logger.error("error %s downloading %s", (exc, self.feed))

    def __unicode__(self):
        return self.title

class BlogPost(models.Model):
    "Represents an entry of a Blog"

    # The blog that generated the entry
    blog = models.ForeignKey(Blog)

    # A unique id for this entry
    uid = models.CharField(max_length=200, default="", null=False)

    # The title of the entry
    title = models.CharField(max_length=200, null=False)

    # The content of the feed
    content = models.TextField(default='', max_length=20000)

    # Santizied HTML
    html = models.TextField(default='')

    # Date related fields.
    creation_date = models.DateTimeField(db_index=True)

    # Date at which the post has been inserted into the database
    insert_date = models.DateTimeField(db_index=True, null=True)

    # Has the entry been published
    published = models.BooleanField(default=False)

    # The link to the entry
    link = models.URLField()

    @property
    def get_title(self):
        return u"BLOG: %s" % self.title

    def get_absolute_url(self):
        return self.link

    def save(self, *args, **kwargs):

        if not self.id:
            # Set the date to current time if missing.
            self.insert_date = self.insert_date or now()

        super(BlogPost, self).save(*args, **kwargs)

    def __unicode__(self):
        return self.title

admin.site.register(Blog)