from django.db import models
from django.conf import settings
import os, logging, datetime
from urllib import request
import feedparser
from django.utils.timezone import utc
from biostar.accounts.models import User
import uuid

logger = logging.getLogger("engine")

MAX_TEXT_LEN = 10000

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
    # Adding field that indicates a remote blog
    remote = models.BooleanField(default=True)

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
            stream = request.urlopen(self.feed)
            text = stream.read().decode("utf-8", errors="replace")
            stream = open(self.fname, 'w', encoding='utf-8')
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

    # Sanitized HTML
    html = models.TextField(default='')

    # Date related fields.
    creation_date = models.DateTimeField(db_index=True)

    # Date at which the post has been inserted into the database
    insert_date = models.DateTimeField(db_index=True, null=True)

    # Has the entry been published
    published = models.BooleanField(default=False)

    # The link to the entry
    link = models.URLField()

    # Posts should be ranked by this.
    rank = models.DateTimeField(db_index=True, null=True)

    @property
    def get_title(self):
        return f"BLOG: {self.title}"

    def get_absolute_url(self):
        return self.link

    def save(self, *args, **kwargs):

        self.insert_date = self.insert_date or now()

        # Set the rank
        self.rank = self.rank or self.insert_date
        #self.html = self.hmtl or
        # SET THE HTML
        #self.html = ''
        self.uid = self.uid or get_uuid(10)

        super(BlogPost, self).save(*args, **kwargs)

    def __str__(self):
        return self.title

