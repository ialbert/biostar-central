from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime
from django import forms
from django.db import models
from django.conf import settings
from django.contrib import admin
from django.contrib.auth import get_user_model
from django.utils.timezone import utc
from taggit.managers import TaggableManager
import reversion


# Obtain the user model
User = get_user_model()

logger = logging.getLogger(__name__)


class Post(models.Model):
    "Represents a post in Biostar"

    tags = TaggableManager()

    # Post statuses.
    OPEN, CLOSED, DELETED = range(3)
    STATUS_CHOICES = [(OPEN, "Open"), (CLOSED, "Closed"), (DELETED, "Deleted")]

    # Question types.
    QUESTION, ANSWER, COMMENT, JOB, FORUM = range(5)
    TYPE_CHOICES = [(QUESTION,"Question"), (ANSWER, "Answer"), (COMMENT, "Comment"), (JOB, "Job"), (FORUM, "Forum")]

    title = models.CharField(max_length=255)

    # The user that originally created the post.
    author = models.ForeignKey(User)

    # The user that edited the post most recently.
    lastedit_user = models.ForeignKey(User, related_name='editor')

    # Indicates the information value of the post.
    rank = models.FloatField(default=0, blank=True)

    # Post status: open, closed, deleted.
    status = models.IntegerField(choices=STATUS_CHOICES, default=OPEN)

    # The type of the post: question, answer, comment.
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)

    # Number of upvotes for the post
    vote_count = models.IntegerField(default=0, blank=True, db_index=True)

    # The number of views for the post.
    view_count = models.IntegerField(default=0, blank=True)

    # The number of replies that a post has.
    reply_count = models.IntegerField(default=0, blank=True)

    # Bookmark count.
    book_count = models.IntegerField(default=0)

    # The total score of the thread (used for top level only)
    thread_score = models.IntegerField(default=0, blank=True, db_index=True)

    # Date related fields.
    creation_date = models.DateTimeField(db_index=True)
    lastedit_date = models.DateTimeField(db_index=True)

    # Stickiness of the post.
    sticky = models.BooleanField(default=False, db_index=True)

    # Indicates whether the post has accepted answer.
    has_accepted = models.BooleanField(default=False, blank=True)

    # This will maintain the ancestor/descendant relationship bewteen posts.
    root = models.ForeignKey('self', related_name="descendants", null=True, blank=True)

    # This will maintain parent/child replationships between posts.
    parent = models.ForeignKey('self', null=True, blank=True, related_name='children')

    # This is the sanitized HTML for display.
    html = models.TextField(default='')

    def save(self, *args, **kwargs):

        if not self.id:
            # This runs only once upon object creation.
            self.creation_date = datetime.datetime.utcnow().replace(tzinfo=utc)
            self.lastedit_date = self.creation_date

        super(Post, self).save(*args, **kwargs)

    def __unicode__(self):
        return "%s: %s (%s)" % (self.get_type_display(), self.title, self.id)

# Posts will have revisions.
reversion.register(Post)

# Revision admin setup.
class PostAdmin(reversion.VersionAdmin):
    list_display = ('title', 'type', 'author')
    fieldsets = (
        (None, {'fields': ('title',)}),
        ('Attributes', {'fields': ('type', 'status', 'sticky',)}),
        ('Content', {'fields': ('tags', 'html', )}),
    )
    search_fields = ('title', 'author__name')

admin.site.register(Post, PostAdmin)


class Vote(models.Model):
    # Post statuses.
    UP, DOWN, BOOKMARK, ACCEPT = range(4)
    TYPE_CHOICES = [(UP, "Up"), (DOWN, "Down"), (BOOKMARK, "Bookmark"), (ACCEPT, "Accept")]

    author = models.ForeignKey(User)
    post = models.ForeignKey(Post, related_name='votes')
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)
    creation_date = models.DateTimeField(db_index=True, auto_now=True)