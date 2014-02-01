from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime, reversion
from django.db import models
from django.conf import settings
from django.contrib import admin
from django.utils.timezone import utc
from taggit.managers import TaggableManager
from django.utils.translation import ugettext_lazy as _

logger = logging.getLogger(__name__)

class PostManager(models.Manager):

    def top_level(self, user):
        "Returns posts based on a user type"
        if user.is_moderator:
            query = self.filter(type__in=Post.TOP_LEVEL)
        else:
            query = self.filter(type__in=Post.TOP_LEVEL, status=Post.OPEN)
        return query

class Post(models.Model):
    "Represents a post in Biostar"

    objects = PostManager()

    tags = TaggableManager()

    # Post statuses.
    PENDING, OPEN, CLOSED, DELETED = range(4)
    STATUS_CHOICES = [(PENDING, "Pending"), (OPEN, "Open"), (CLOSED, "Closed"), (DELETED, "Deleted")]

    # Question types.
    QUESTION, ANSWER, COMMENT, JOB, FORUM, PAGE = range(6)
    TYPE_CHOICES = [
        (QUESTION,"Question"), (ANSWER, "Answer"), (COMMENT, "Comment"),
        (JOB, "Job"), (FORUM, "Forum"), (PAGE, "Page"),
    ]

    TOP_LEVEL = set((QUESTION, JOB, FORUM, PAGE))

    title = models.CharField(max_length=255, null=False)

    # The user that originally created the post.
    author = models.ForeignKey(settings.AUTH_USER_MODEL)

    # The user that edited the post most recently.
    lastedit_user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='editor')

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

            if not (self.parent or self.root):
                # Neither of root nor parent are set.
                # This will be a top level post.
                self.root = self.parent = self
                self.type = self.type if self.type in Post.TOP_LEVEL else self.FORUM
            elif self.parent:
                # If the parent is set root must follow the parent.
                self.root = self.parent.root
                self.type = self.ANSWER if self.parent.type in Post.TOP_LEVEL else self.COMMENT
            elif self.root:
                # If only the root is set the parent has follow the root.
                # Note: In general the root should not be directly set.
                # Let the save method figure out the correct root.
                self.parent = self.root
                self.type = self.ANSWER
            else:
                raise Exception("This shouldn't ever happen really")

            self.title = self.parent.title if self.parent else self.title
            self.lastedit_user = self.author
            self.status = self.status or Post.PENDING
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

    author = models.ForeignKey(settings.AUTH_USER_MODEL)
    post = models.ForeignKey(Post, related_name='votes')
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)
    date = models.DateTimeField(db_index=True, auto_now=True)

class SubscriptionManager(models.Manager):

    def get_subs(self, post):
        "Returns all suscriptions for a post"
        return self.filter(post=post.root).select_related("user")

# This contains the notification types.
from biostar.const import LOCAL_MESSAGE, MESSAGING_TYPE_CHOICES

class Subscription(models.Model):
    "Connects a post to a user"

    class Meta:
        unique_together = (("user", "post"),)

    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name=_("User"),  db_index=True)
    post = models.ForeignKey(Post, verbose_name=_("Post"), related_name="subs", db_index=True)
    type = models.IntegerField(choices=MESSAGING_TYPE_CHOICES, default=LOCAL_MESSAGE, db_index=True)
    date = models.DateTimeField(_("Date"), db_index=True)

    objects = SubscriptionManager()

    def __unicode__(self):
        return "%s to %s" % (self.user.name, self.post.title)

    @staticmethod
    def create(post, user):
        "Creates a subscription of a user to a post"
        root = post.root
        if Subscription.objects.filter(post=root, user=user).count() == 0:
            sub = Subscription(post=root, user=user)
            sub.date = datetime.datetime.utcnow().replace(tzinfo=utc)
            sub.save()


# Admin interface for subscriptions
class SubscriptionAdmin(admin.ModelAdmin):
    search_fields = ('user__name', 'user__email')
    list_select_related = ["user", "post"]

admin.site.register(Subscription, SubscriptionAdmin)

