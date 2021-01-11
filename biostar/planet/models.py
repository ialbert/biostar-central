from django.db import models
from django.conf import settings
import os, logging, datetime
from urllib import request
import feedparser
from django.utils.timezone import utc
import uuid

logger = logging.getLogger(__name__)


def now():
    return datetime.datetime.utcnow().replace(tzinfo=utc)


def abspath(*args):
    """
    Generates absolute paths
    """
    return os.path.abspath(os.path.join(*args))


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


class Blog(models.Model):
    """
    Represents a blog
    """
    title = models.CharField(max_length=255, default="")
    desc = models.TextField(default='', blank=True)
    feed = models.URLField()
    link = models.URLField()
    active = models.BooleanField(default=True)
    list_order = models.IntegerField(default=0)

    @property
    def fname(self):
        fname = abspath(settings.PLANET_DIR, f"{self.id}.xml")
        return fname

    def parse(self):
        try:
            doc = feedparser.parse(self.feed)
        except Exception as exc:
            logger.error(f"Error parsing feed. {exc}")
            doc = None
        return doc

    def download(self):
        try:
            text = request.urlopen(self.feed).read().decode()
            stream = open(self.fname, 'wt')
            stream.write(text)
            stream.close()
        except Exception as exc:
            logger.error(f"Error downloading {exc}")

    def __str__(self):
        return self.title


class BlogPost(models.Model):
    "Represents an entry of a Blog"

    # The blog that generated the entry
    blog = models.ForeignKey(Blog, on_delete=models.CASCADE)

    # A unique id for this entry
    uid = models.CharField(max_length=200, unique=True)

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
        return f"BLOG: {self.title}"

    def get_absolute_url(self):
        return self.link

    def save(self, *args, **kwargs):

        self.insert_date = self.insert_date or now()
        self.uid = self.uid or get_uuid(10)
        super(BlogPost, self).save(*args, **kwargs)

    def __str__(self):
        return self.title

