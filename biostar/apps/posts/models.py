from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime
from django import forms
from django.db import models
from django.conf import settings
from django.contrib import admin
from django.contrib.auth import get_user_model
from django.utils.timezone import utc

# Obtain the user model
User = get_user_model()

logger = logging.getLogger(__name__)


class Post(models.Model):

    # Post statuses.
    OPEN, CLOSED, DELETED = range(3)
    STATUS_CHOICES = [(OPEN, "Open"), (CLOSED, "Closed"), (DELETED, "Deleted")]

    # Question types.
    QUESTION, ANSWER, COMMENT, JOB, FORUM = range(5)
    TYPE_CHOICES = [(QUESTION,"Question"), (ANSWER, "Answer"), (COMMENT, "Comment"), (JOB, "Job"), (FORUM, "Forum")]

    title = models.TextField(max_length=200)

    # The tag string is the canonical form of the tags
    tag_string = models.CharField(max_length=200, blank=True,)

    # Post attributes
    views = models.IntegerField(default=0, blank=True, db_index=True)
    score = models.IntegerField(default=0, blank=True, db_index=True)
    full_score = models.IntegerField(default=0, blank=True, db_index=True)

    creation_date = models.DateTimeField(db_index=True)
    lastedit_date = models.DateTimeField(db_index=True)

    # The user that originally created the post.
    author = models.ForeignKey(User)

    # The user that edited the post most recently.
    lastedit_user = models.ForeignKey(User, related_name='editor')

    # The number of replies that a post has.
    reply_count = models.IntegerField(default=0, blank=True)

    # Bookmark count.
    book_count = models.IntegerField(default=0, blank=True)

    # Stickiness of the post.
    sticky = models.IntegerField(default=0, db_index=True)

    # Indicates whether the post has accepted answer.
    has_accepted = models.BooleanField(default=False, blank=True)

    # Post status: open, closed, deleted.
    status = models.IntegerField(choices=STATUS_CHOICES, default=OPEN)

    # The type of the post: question, answer, comment.
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)

    # Some posts have links to other content.
    link_out = models.URLField(default='', blank=True)

    # This will maintain the ancestor/descendant relationship bewteen posts.
    root = models.ForeignKey('self', related_name="descendants", null=True, blank=True)

    # This will maintain parent/child replationships between posts.
    parent = models.ForeignKey('self', null=True, blank=True, related_name='children')

    def save(self, *args, **kwargs):

        if not self.id:
            # This runs only once upon object creation.
            self.creation_date = datetime.datetime.utcnow().replace(tzinfo=utc)
            self.lastedit_date = self.creation_date

        super(Post, self).save(*args, **kwargs)

    def __unicode__(self):
        return "%s: %s (%s)" % (self.get_type_display(), self.title, self.id)

class PostBody(models.Model):
    "The post stores each revision"
    post = models.OneToOneField(Post, related_name="body")
    author = models.ForeignKey(User)

    # Title and tag_strings are duplicated because they also serve as revision.
    title = models.TextField(max_length=200)

    # The underlying Markdown
    content = models.TextField(default='', max_length=25000)

    # This is the sanitized HTML for display.
    html = models.TextField(default='')

